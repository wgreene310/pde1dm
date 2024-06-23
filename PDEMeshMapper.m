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
%   Copyright (C) 2016-2024 William H. Greene
classdef PDEMeshMapper < handle
  
  %% Map a solution and derivatives from FEM nodes to arbitrary points
  %% in the region.
  
  properties
    destNodeIndex
  end
  
  methods
    
    function obj = PDEMeshMapper(srcMesh, destMesh, numElemNodes, sF, dSf)
      obj.numElemNodes = numElemNodes;
      obj.sF=sF;
      obj.dSf = dSf;

      obj.destMesh = destMesh;
      numDest = length(destMesh);
      obj.destMeshFirstLastNodeIndex = zeros(numDest,2);
      obj.destMeshParamVals = zeros(numDest,1);
      obj.N = zeros(numElemNodes, numDest);
      obj.dNdx = zeros(numElemNodes, numDest);
      
      if 1
        obj.setSrcMesh(srcMesh);
      else
        obj.srcMesh = srcMesh;
        numSrc=length(srcMesh);
        nen = obj.numElemNodes;
        firstNode = @(elemIndex) (nen-1)*(elemIndex-1) + 1;
        getElemIndex=@(nodeId) nodeId/(nen-1);
        numDestMeshNodes = 0;
        for i=1:numDest
          xd=destMesh(i);
          id=find(srcMesh<=xd, 1, 'last');
          if isempty(id)
            indLeft=1;
            indRight = nen;
          elseif(id==numSrc)
            indLeft=id-nen+1;
            indRight = numSrc;
          else
            ei=ceil(getElemIndex(id));
            indLeft = firstNode(ei);
            indRight=indLeft+nen-1;
          end
          obj.destMeshFirstLastNodeIndex(i,:) = [indLeft indRight];
          obj.destNodeIndex = [obj.destNodeIndex indLeft:indRight];
          
          numDestMeshNodes = numDestMeshNodes+numElemNodes;
          s=obj.calcXi(srcMesh(indLeft), srcMesh(indRight), xd);
          obj.N(:,i)=obj.sF(s);
          ei=obj.destMeshFirstLastNodeIndex(i,:);
          jac = 2. / (srcMesh(ei(2))-srcMesh(ei(1)));
          obj.dNdx(:,i) = jac*obj.dSf(s);
          obj.destMeshParamVals(i) = s;
        end
        obj.allocTmp = true;
        obj.destNodeIndex = unique(obj.destNodeIndex);
      end
    end
    
    function setSrcMesh(self, srcMesh)
      self.srcMesh = srcMesh;
      numSrc=length(srcMesh);
      nen = self.numElemNodes;
      firstNode = @(elemIndex) (nen-1)*(elemIndex-1) + 1;
      getElemIndex=@(nodeId) nodeId/(nen-1);
      numDestMeshNodes = 0;
      for i=1:length(self.destMesh)
        xd=self.destMesh(i);
        id=find(srcMesh<=xd, 1, 'last');
        if isempty(id)
          indLeft=1;
          indRight = nen;
        elseif(id==numSrc)
          indLeft=id-nen+1;
          indRight = numSrc;
        else
          ei=ceil(getElemIndex(id));
          indLeft = firstNode(ei);
          indRight=indLeft+nen-1;
        end
        self.destMeshFirstLastNodeIndex(i,:) = [indLeft indRight];
        self.destNodeIndex = [self.destNodeIndex indLeft:indRight];
        
        numDestMeshNodes = numDestMeshNodes+nen;
        s=self.calcXi(srcMesh(indLeft), srcMesh(indRight), xd);
        self.N(:,i)=self.sF(s);
        ei=self.destMeshFirstLastNodeIndex(i,:);
        jac = 2. / (srcMesh(ei(2))-srcMesh(ei(1)));
        self.dNdx(:,i) = jac*self.dSf(s);
        self.destMeshParamVals(i) = s;
      end
      self.allocTmp = true;
      self.destNodeIndex = unique(self.destNodeIndex);
    end
  
    function destU=mapFunction(self,srcU)
      destU=self.mapFunctionImpl(srcU, false);
    end
    
    function destU=mapFunctionDer(self,srcU)
      destU=self.mapFunctionImpl(srcU, true);
    end
    
    function uOut=mapFunctionImpl(self,srcU, calcDeriv)
      numDepVars = size(srcU,1);
      numDest = length(self.destMesh);
      if self.allocTmp
        self.destU = zeros(numDepVars, numDest, 'like', srcU);
        self.allocTmp=false;
      end
      for i=1:numDest
        ei=self.destMeshFirstLastNodeIndex(i,:);
        if(calcDeriv)
          self.destU(:,i) = srcU(:,ei(1):ei(2))*self.dNdx(:,i);
        else
          self.destU(:,i) = srcU(:,ei(1):ei(2))*self.N(:,i);
        end
      end
      uOut=self.destU;
    end
    
  end % methods
  
  methods(Static)
    function xid=calcXi(xim1, xi, x)
      xid= 2 * (x - xim1) / (xi - xim1) - 1;
    end
  end % methods
  
  properties(Access=private)
    srcMesh, destMesh, numElemNodes;
    destMeshFirstLastNodeIndex, destMeshParamVals;
    sF, dSf, N, dNdx;
    allocTmp, destU; % temporary for use in mapFunctionImpl
  end
  
end

