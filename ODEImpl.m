classdef ODEImpl
  
  properties
    numODEDOFs;
    y0Ode;
  end
  
  methods
    
    function obj = ODEImpl(xmesh, odeFunc, odeICFunc, odeMesh)
      obj.odeFunc = odeFunc;
      obj.odeICFunc = odeICFunc;
      obj.odeMesh = odeMesh;
      obj.meshMapper = PDEMeshMapper(xmesh, odeMesh);
      obj.y0Ode = odeICFunc();
      obj.numODEVariables = length(obj.y0Ode);
      obj.numODEEquations = 0; % FIXME
      obj.numODEDOFs = obj.numODEVariables + obj.numODEEquations;
    end
    
    function f=calcODEResidual(self, time, u2, up2, f2, v, vDot)
      mm = self.meshMapper;
      uOde=mm.mapFunction(u2);
      upOde=mm.mapFunction(up2);
      dudxOde=mm.mapFunctionDer(u2);
      fluxOde = mm.mapFunction(f2);
      dupdxOde=mm.mapFunctionDer(up2);
      f = callVarargFunc(self.odeFunc, ...
        {time, v, vDot, self.odeMesh, uOde, dudxOde, fluxOde, upOde, dupdxOde});
    end
    
   function [dFdv,dFdvDot]=calcDOdeDv(self, time, u2, up2, f2, v, vDot)
      mOde = self.numODEVariables;
      mm = self.meshMapper;
      uOde=mm.mapFunction(u2);
      upOde=mm.mapFunction(up2);
      dudxOde=mm.mapFunctionDer(u2);
      fluxOde = mm.mapFunction(f2);
      dupdxOde=mm.mapFunctionDer(up2);
      f0 = callVarargFunc(self.odeFunc, ...
        {time, v, vDot, self.odeMesh, uOde, dudxOde, fluxOde, upOde, dupdxOde});
      if size(f0,1) ~= mOde
        error('Number of rows returned from ODE function must be %d.\n', ...
          mOde);
      end
      sqrtEps = sqrt(eps);
      dFdv = zeros(mOde, mOde);
      for i=1:mOde
        vSave = v(i);
        h = sqrtEps * max(vSave, 1);
        v(i) = v(i) + h;
        fp = callVarargFunc(self.odeFunc, ...
          {time, v, vDot, self.odeMesh, uOde, dudxOde, fluxOde, upOde, dupdxOde});
        dFdv(:,i) = (fp-f0)/h;
        v(i) = vSave;
      end
      dFdvDot = zeros(mOde, mOde);
      for i=1:mOde
        vDotSave = vDot(i);
        h = sqrtEps * max(vDotSave, 1);
        vDot(i) = vDot(i) + h;
        fp = callVarargFunc(self.odeFunc, ...
          {time, v, vDot, self.odeMesh, uOde, dudxOde, fluxOde, upOde, dupdxOde});
        dFdvDot(:,i) = (fp-f0)/h;
        vDot(i) = vDotSave;
      end
   end
 
    function dFdu=calcDOdeDu(self, time, u2, up2, f2, v, vDot)
      [numDepVars,numNodes] = size(u2);
      numFEMEqns = numDepVars*numNodes;
      mm = self.meshMapper;
      uOde=mm.mapFunction(u2);
      upOde=mm.mapFunction(up2);
      dudxOde=mm.mapFunctionDer(u2);
      fluxOde = mm.mapFunction(f2);
      dupdxOde=mm.mapFunctionDer(up2);
      f0 = callVarargFunc(self.odeFunc, ...
        {time, v, vDot, self.odeMesh, uOde, dudxOde, fluxOde, upOde, dupdxOde});
      if size(f0,1) ~= self.numODEVariables
        error('Number of rows returned from ODE function must be %d.\n', ...
          self.numODEVariables);
      end
      sqrtEps = sqrt(eps);
      dFdu = zerosLike(self.numODEVariables, numDepVars, numNodes, u2);
      for i=1:numNodes
        for j=1:numDepVars
          usave = u2(j,i);
          h = sqrtEps * max(usave, 1);
          u2(j,i) = u2(j,i) + h;
          uOde=mm.mapFunction(u2);
          dudxOde=mm.mapFunctionDer(u2);
          f = callVarargFunc(self.odeFunc, ...
            {time, v, vDot, self.odeMesh, uOde, dudxOde, fluxOde, upOde, dupdxOde});
          dFdu(:,j,i) = (f-f0)/h;
          u2(j,i) = usave;
        end
      end
      dFdu = reshape(dFdu, self.numODEVariables, numFEMEqns);
    end    

  end
  
  methods(Access=private)
  end
  
  properties(Access=private)
    odeFunc, odeICFunc, odeMesh;
    numODEVariables, numODEEquations;
    meshMapper;
  end
  
end

