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
% Copyright (C) 2016-2021 William H. Greene

function [sol,varargout] = pde1dM (m, pde,ic,bc,xmesh,t,varargin)

if(~ismember(nargin,[6 7 9 10]))
  error('pde1d:nrhs', ...
  'Illegal number of arguments passed to pde1dM');
end
hasODE = nargin > 7;

p=inputParser;
validOnOff = @(x) any(validatestring(x,{'on','off'}));
p.addRequired('m', @validM);
validHandle = @(x) isa(x, 'function_handle');
p.addRequired('pdeFunc', validHandle);
p.addRequired('icFunc', validHandle);
p.addRequired('bcFunc', validHandle);
p.addRequired('xmesh', @validX);
p.addRequired('t', @validT);
if hasODE
  p.addOptional('odeFunc', validHandle);
  p.addOptional('odeICFunc', validHandle);
  validOdeMesh = @(x) isreal(x);
  p.addOptional('odeMesh', validOdeMesh);
end
p.addParameter('RelTol', 1e-3);
p.addParameter('AbsTol', 1e-6);
p.addParameter('Stats', 'off',validOnOff);
p.addParameter('BDF', 'off');
p.addParameter('MaxOrder', 5);
p.addParameter('MaxStep', []);
p.addParameter('Vectorized', 'off',validOnOff);
p.addParameter('EqnDiagnostics', 0);
p.addParameter('ICDiagnostics', 0);
p.addParameter('PolyOrder', 1);
p.addParameter('ViewMesh', 1);
p.addParameter('DiagonalMassMatrix', 'off', validOnOff);
validPts = @(n) isscalar(n) && n>0;
p.addParameter('NumIntegrationPoints', 2, validPts);
p.addParameter('AnalyticalJacobian', 'off', validOnOff);
p.addParameter('InitialSlope', [], validHandle);
validInitFunc = @(x) isempty(x) || validHandle(x);
p.addParameter('eqnDiagnosticsInitFunc', [], validInitFunc);
p.parse(m,pde,ic,bc,xmesh,t,varargin{:});
%p.Results

pdeOpts = PDEOptions;
pdeOpts.hasODE = hasODE;
pdeOpts.icDiagnostics=p.Results.ICDiagnostics;
eqnDiagnostics=p.Results.EqnDiagnostics;
pdeOpts.eqnDiagnostics=eqnDiagnostics;
pdeOpts.useDiagMassMat = strcmpi(p.Results.DiagonalMassMatrix, 'on');
pdeOpts.vectorized = strcmpi(p.Results.Vectorized, 'on');
pdeOpts.numIntegrationPoints = p.Results.NumIntegrationPoints;
pdeOpts.analyticalJacobian = strcmpi(p.Results.AnalyticalJacobian, 'on');
pdeOpts.initialSlope = p.Results.InitialSlope;
pdeOpts.eqnDiagnosticsInitFunc = p.Results.eqnDiagnosticsInitFunc;

if  hasODE
  pdeImpl = PDE1dImpl(m, pde, ic, bc, xmesh, t, pdeOpts, ...
  p.Results.odeFunc, p.Results.odeICFunc, ...
  p.Results.odeMesh);
else
  pdeImpl = PDE1dImpl(m, pde, ic, bc, xmesh, t, pdeOpts);
end


if(~eqnDiagnostics)
  odeOpts=odeset('RelTol', p.Results.RelTol, ...
    'AbsTol', p.Results.AbsTol, ...
    'Stats', p.Results.Stats, ...
    'BDF', p.Results.BDF, ...
    'MaxOrder', p.Results.MaxOrder, ...
    'MaxStep', p.Results.MaxStep);
  if(hasODE)
    [u,uODE] = pdeImpl.solveTransient(odeOpts);
    varargout{1} = uODE;
  else
    u = pdeImpl.solveTransient(odeOpts);
  end
  npde=pdeImpl.numDepVars;
  if(npde==1)
    sol = u;
  else
    sol = toSol3(u, npde, length(xmesh));
  end
else
  pdeImpl.testFuncs;
  sol=[];
  varargout{1} = [];
end
end

function sol=toSol3(u, npde, nx)
nt=size(u, 1);
sol = zeros(nt, nx, npde);
for i = 1:npde
  sol(:,:,i) = u(:,i:npde:end);
end
end

function valid=validM(m)
valid = isscalar(m) && isreal(m) && (m==0 || m==1 || m==2);
if(~valid)
  error('pde1d:invalid_m_val', 'Illegal value of m');
end
end

function valid=validX(x)
valid =  isreal(x) && length(x)>=2;
if(~valid)
  error('pde1d:mesh_type', 'Illegal value of xmesh');
end
end

function valid=validT(x)
valid =  isreal(x) && length(x)>=2;
if(~valid)
  error('pde1d:time_type', 'Illegal value of tspan');
end
end
