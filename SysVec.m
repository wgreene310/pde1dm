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
% Copyright (C) 2020-2021 William H. Greene

classdef SysVec
  
  %% Represent a solution vector with finite element
  %% and ODE DOFs.
  
  methods
    function obj = SysVec(v, numFEMNodes, numDepVars)
      obj.dataVector = v;
      obj.numFEMNodes=numFEMNodes;
      if nargin > 2
        obj.numDepVars=numDepVars;
      else
        obj.numDepVars=1;
      end
      obj.numFEMEqns=obj.numFEMNodes*obj.numDepVars;
    end
    
    function v=new(self, data)
      v=SysVec(data, self.numFEMNodes, self.numDepVars);
    end
    
    function v=fem(self)
      v = self.dataVector(1:self.numFEMEqns);
    end
    
    function v2=fem2(self)
      if self.numDepVars==1
        v2=self.dataVector(1:self.numFEMEqns)';
      else
        v2 = reshape(self.dataVector(1:self.numFEMEqns), self.numDepVars, ...
          self.numFEMNodes);
      end
    end
    
    function vode=ode(self)
      vode = self.dataVector(self.numFEMEqns+1:end);
    end
    
    function v=all(self)
      v=self.dataVector;
    end
  end
  
  properties(Access=private)
    dataVector, numFEMNodes, numDepVars;
    numFEMEqns;
  end
end

