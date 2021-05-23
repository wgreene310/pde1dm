%% provide replacement for pdepe function in octave
function [sol,varargout] = pdepe(m, pde,ic,bc,xmesh,t,varargin)
[sol,varargout{1:nargout-1}] = pde1dM (m, pde,ic,bc,xmesh,t,varargin{:});
end