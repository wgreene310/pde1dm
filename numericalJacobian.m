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
% Copyright (C) 2020 William H. Greene

function J=numericalJacobian(f, S, y, t, groupList)
      J=S;
      [m,n]=size(S);
      r0=f(t,y);
      sqrtEps = sqrt(eps);
      ySave = zeros(n,1);
      numGroups=max(groupList);
      %numGroups=1; % TESTING !!!
      for g=1:numGroups
        %gCols = find(groupList,g);
        gCols = find(groupList==g);
        %prtShortVec(gCols, 'gCols');
        ySave(gCols) = y(gCols);
        h = sqrtEps*max(ySave(gCols),1);
        y(gCols) = y(gCols) + h;
        %prtShortVec(y, 'y');
        dr=f(t,y)-r0;
        %dr(gCols) = dr(gCols)./h;
        for i=1:length(gCols)
          c=gCols(i);
          rowsI=find(S(:,c));
          J(rowsI,c) = dr(rowsI)/h(i);
        end
        y(gCols) = ySave(gCols);
      end
    end