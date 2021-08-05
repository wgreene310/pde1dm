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

classdef PDEOptions
  
  %% All the user-settable options for the PDE solver.
  
  properties
    icDiagnostics, eqnDiagnostics, useDiagMassMat;
    vectorized, numIntegrationPoints, hasODE, analyticalJacobian;
    initialSlope, eqnDiagnosticsInitFunc, isOctave;
    useInternalNumJac;
    testFunctionDOFMap;
    polyOrder;
  end
  
  methods
    function obj = PDEOptions()
      obj.icDiagnostics=0;
      obj.eqnDiagnostics=0;
      obj.eqnDiagnosticsInitFunc=[];
      obj.useDiagMassMat = false;
      obj.vectorized = false;
      obj.numIntegrationPoints = [];
      obj.hasODE = false;
      obj.analyticalJacobian = false;
      obj.initialSlope = [];
      obj.isOctave = exist('OCTAVE_VERSION', 'builtin');
      if obj.isOctave
        obj.useInternalNumJac=true;
      else
        obj.useInternalNumJac=false;
      end
      obj.testFunctionDOFMap = [];
      obj.polyOrder=1;
    end
  end
end

