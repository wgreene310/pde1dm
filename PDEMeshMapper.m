%  
%   This library is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 3
%   of the License, or (at your option) any later version.
%  
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   General Public License for more details.
%  
%   You should have received a copy of the GNU General Public License
%   along with this library; if not, visit
%   http://www.gnu.org/licenses/gpl.html or write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
% 
%   Copyright (C) 2016-2021 William H. Greene
classdef PDEMeshMapper
  
  %% Map a solution and derivatives from FEM nodes to arbitrary points
  %% in the region.
  
  properties
    destNodeIndex
  end
  
  methods
    
    function obj = PDEMeshMapper(srcMesh, destMesh, numElemNodes, sF, dSf)
      obj.srcMesh = srcMesh;
      numSrc=length(srcMesh);
      obj.destMesh = destMesh;
      numDest = length(destMesh);
      obj.destMeshFirstLastNodeIndex = zeros(numDest,2);
      obj.destMeshParamVals = zeros(numDest,1);
      obj.N = zeros(numElemNodes, numDest);
      obj.dNdx = zeros(numElemNodes, numDest);
      numDestMeshNodes = 0;
      for i=1:numDest
        xd=destMesh(i);
        id=find(srcMesh<=xd, 1, 'last');
        if(id==1)
          indLeft=1;
          indRight = numElemNodes;
        elseif(id==numSrc)
          indLeft=id-numElemNodes+1;
          indRight = numSrc;
        else
          indLeft = id;
          indRight=id+numElemNodes-1;
        end
        obj.destMeshFirstLastNodeIndex(i,:) = [indLeft indRight];
        obj.destNodeIndex = [obj.destNodeIndex indLeft:indRight];
        
        numDestMeshNodes = numDestMeshNodes+numElemNodes;
        s=obj.calcXi(srcMesh(indLeft), srcMesh(indRight), xd);
        obj.N(:,i)=sF(s);
        ei=obj.destMeshFirstLastNodeIndex(i,:);
        jac = 2. / (srcMesh(ei(2))-srcMesh(ei(1)));
        obj.dNdx(:,i) = jac*dSf(s);
        obj.destMeshParamVals(i) = s;
      end
      obj.destNodeIndex = unique(obj.destNodeIndex);
    end
  
    function destU=mapFunction(self,srcU)
      destU=self.mapFunctionImpl(srcU, false);
    end
    
    function destU=mapFunctionDer(self,srcU)
      destU=self.mapFunctionImpl(srcU, true);
    end
    
    function destU=mapFunctionImpl(self,srcU, calcDeriv)
      numDepVars = size(srcU,1);
      numDest = length(self.destMesh);
      destU = zerosLike(numDepVars, numDest, srcU);
      for i=1:numDest
        ei=self.destMeshFirstLastNodeIndex(i,:);
        if(calcDeriv)
          destU(:,i) = srcU(:,ei(1):ei(2))*self.dNdx(:,i);
        else
          destU(:,i) = srcU(:,ei(1):ei(2))*self.N(:,i);
        end       
      end
    end
    
  end % methods
  
  methods(Static)
    function xid=calcXi(xim1, xi, x)
      xid= 2 * (x - xim1) / (xi - xim1) - 1;
    end
  end % methods
  
  properties(Access=private)
    srcMesh, destMesh;
    destMeshFirstLastNodeIndex, destMeshParamVals;
    N, dNdx;
  end
  
end

