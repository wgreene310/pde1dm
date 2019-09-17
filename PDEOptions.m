classdef PDEOptions
  properties
    icDiagnostics, eqnDiagnostics, addLagMultVector, useDiagMassMat;
    vectorized, numIntegrationPoints, hasODE, analyticalJacobian;
  end
  
  methods
    function obj = PDEOptions()
      obj.icDiagnostics=0;
      obj.eqnDiagnostics=0;
      obj.addLagMultVector = false;
      obj.useDiagMassMat = false;
      obj.vectorized = false;
      obj.numIntegrationPoints = 2;
      obj.hasODE = false;
      obj.analyticalJacobian = false;
    end
  end
end

