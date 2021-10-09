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

function xmesh=pdeRefineMesh(xmesh, polyOrder)
if polyOrder<2
  return;
end
% add intermediate nodes for higher-order elems
numNodesOrig = length(xmesh);
numElems=numNodesOrig-1;
elemLen=xmesh(2:end)-xmesh(1:end-1);
del=elemLen/polyOrder;
numNodesNew=numElems*polyOrder+1;
xn=zeros(1,numNodesNew);
ii=1:polyOrder:numNodesNew-polyOrder;
xn([ii end])=xmesh;
for i=2:polyOrder
  ii=ii+1;
  xn(ii)=xmesh(1:end-1)+(i-1)*del;
end
xmesh=xn;
end