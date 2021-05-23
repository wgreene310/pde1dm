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
%   Copyright (C) 2016-2020 William H. Greene
classdef PDEMeshMapper
  
  %% Map a solution and derivatives from FEM nodes to arbitrary points
  %% in the region.
  
  properties
    srcMesh, destMesh;
    destMeshElemIndex, destMeshParamVals;
  end
  
  methods
    
    function obj = PDEMeshMapper(srcMesh, destMesh)
      obj.srcMesh = srcMesh;
      numSrc=length(srcMesh);
      obj.destMesh = destMesh;
      numDest = length(destMesh);
      obj.destMeshElemIndex = zeros(numDest,1);
      obj.destMeshParamVals = zeros(numDest,1);
      for i=1:numDest     
        xd=destMesh(i);
        id=find(srcMesh<=xd, 1, 'last');
        if(id==1)
          indRight = 2;
        elseif(id==numSrc)
          indRight = numSrc;
        else
          indRight=id+1;
        end
        obj.destMeshElemIndex(i) = indRight;
        s=obj.calcXi(srcMesh(indRight-1), srcMesh(indRight), xd);
        obj.destMeshParamVals(i) = s;
      end
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
        s = self.destMeshParamVals(i);
        ei=self.destMeshElemIndex(i);
        if(calcDeriv)
          dN = self.dShapeLine2(s);
          jac = 2. / (self.srcMesh(ei)-self.srcMesh(ei-1));
          destU(:,i) = jac*(dN(1)*srcU(:,ei-1) + dN(2)*srcU(:,ei));
        else
          N = self.shapeLine2(s);
          destU(:,i) = N(1)*srcU(:,ei-1) + N(2)*srcU(:,ei);
        end
      
      end
    end
  
  end % methods
  
  methods(Static)
    function xid=calcXi(xim1, xi, x)
      xid= 2 * (x - xim1) / (xi - xim1) - 1;
    end
    function shp = shapeLine2( r )
      shp = [(1-r)/2 (1+r)/2]';
    end
    function ds = dShapeLine2( r )
      ds = [-.5 .5]'*ones(1,length(r));
    end
  end % methods
  
end

