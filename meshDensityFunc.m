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
% Copyright (C) 2022-2023 William H. Greene

function rho = meshDensityFunc(x,u, varargin)
npde=size(u,1);
nx=length(x);

p=inputParser;
p.addRequired('x', ...
  @(n) validateattributes(n, {'numeric'}, {}));
p.addRequired('u', ...
  @(n) validateattributes(n, {'numeric'}, {}));
p.addOptional('monitor', 3, ...
  @(n) validateattributes(n, {'numeric'}, {}));
p.addOptional('debug', 0, ...
  @(n) validateattributes(n, {'numeric'}, {}));
p.parse(x, u, varargin{:});
monitor=p.Results.monitor;
debug=p.Results.debug;

if debug
  pruneTol=1e-14;
  prtForm='%.10g ';
  prtShortVec(x, 'meshDensityFunc:x', 1, prtForm, pruneTol);
  prtMat(u, 'meshDensityFunc:v', 1, prtForm, pruneTol);
end

duDxx = zeros(nx,1);
uTmp = zeros(nx,1);

% huang eqn 1.22 for second derivative
for i=1:npde
    jj=2:(nx-1);
    % central diff for interior points
    uTmp(jj)=2./(x(jj+1)-x(jj-1)).*(...
      (u(i,jj+1)-u(i,jj))'./(x(jj+1)-x(jj)) ...
      -(u(i,jj)-u(i,jj-1))'./(x(jj)-x(jj-1)));
  uTmp(1) = 2*((x(2)-x(1))*(u(i,3)-u(i,1))-(x(3)-x(1))*(u(i,2)-u(i,1))) ...
    /((x(3)-x(1))*(x(2)-x(1))*(x(3)-x(2)));
  uTmp(nx) = 2*((x(nx-1)-x(nx))*(u(i,nx-2)-u(i,nx)) ...
    -(x(nx-2)-x(nx))*(u(i,nx-1)-u(i,nx))) ...
    /((x(nx-2)-x(nx))*(x(nx-1)-x(nx))*(x(nx-2)-x(nx-1)));
  duDxx = duDxx + uTmp.^2;
end

if debug>2
  figure; plot(x,uTmp); title v2; grid;
end



if (monitor==4)
  % compute the arclength mesh density function
  j=1:nx;
  rho=sqrt(1+duDxx(j));
elseif (monitor==5)
  % compute the curvature mesh density function
  j=1:nx;
  rho=(1+duDxx(j))^0.25;
else  
  % compute alpha
  k=1; m=1;
  gamma = 1.0/(2*(k-m+1)+1);
  alpha = eps;
  jj=2:nx;
  alpha = alpha + 0.5*sum((duDxx(jj).^gamma+duDxx(jj-1).^gamma).*(x(jj)-x(jj-1)));
  alpha = (alpha/(x(nx)-x(1)))^(1/gamma);
  
  rho = (1+(1/alpha)*duDxx).^gamma;
  
end

if debug
  prtShortVec(rho, 'meshDensityFunc:rho', 1, prtForm, pruneTol);
end


end 

