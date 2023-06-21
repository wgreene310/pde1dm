%% provide replacement for pdepe function in octave
function [sol,varargout] = pdepe(m, pde,ic,bc,xmesh,t,opts)
if nargin < 6
  error('pdepe requires at least six input arguments');
end
if nargin < 7
  opts=struct;
end
opts.automaticBCAtCenter=true;
[sol,varargout{1:nargout-1}] = pde1dm (m, pde,ic,bc,xmesh,t,opts);
end