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
% Copyright (C) 2016-2023 William H. Greene
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

function resX=mmPDEResidual(monitor,mmpde,tau,x,xt,rho)

nx=length(x);
resX=zeros(nx,1);

if monitor==0
  % keep mesh fixed
  resX = xt(:);
else
  % moving mesh
  dXi = nx-1; % xi ranges from zero to one
  j=2:(nx-1); % interior nodes
  if mmpde==4
    resX(j) = (tau*(xt(j+1)-xt(j))+(x(j+1)-x(j))).*(rho(j+1)+rho(j)) ...
             -(tau*(xt(j)-xt(j-1))+(x(j)-x(j-1))).*(rho(j)+rho(j-1));    
  elseif mmpde==6  
    resX(j) = (tau*(xt(j+1)-xt(j))+(x(j+1)-x(j)).*(rho(j+1)+rho(j))) ...
             -(tau*(xt(j)-xt(j-1))+(x(j)-x(j-1)).*(rho(j)+rho(j-1))); 
  elseif mmpde==5   
    resX(j) = tau*xt(j)-0.5*dXi*dXi*((x(j+1)-x(j)).*(rho(j+1)+rho(j)) ...
                                    -(x(j)-x(j-1)).*(rho(j)+rho(j-1)));   
  else % for modified MMPDE5 (default)    
    resX(j) = tau*xt(j)-0.5*dXi*dXi/rho(j)*((x(j+1)-x(j)).*(rho(j+1)+rho(j)) ...
                                           -(x(j)-x(j-1)).*(rho(j)+rho(j-1)));   
  end
  
  % boundary points don't move  
  resX(1) = xt(1)-0.0;
  resX(nx) = xt(nx)-0.0;
end

end