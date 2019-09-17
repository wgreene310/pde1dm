classdef PDEMeshMapper
  %PDEMeshMapper 
  
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

