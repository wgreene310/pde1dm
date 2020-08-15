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

function groupList=findJacobianGroups(a)
      [m,n]=size(a);
      groupList=zeros(n,1);
      currentGroup=1;
      groupEqns=zeros(m,1);
      while true
        colsNotInGroup=find(groupList==0);
        if isempty(colsNotInGroup)
          break;
        end
        i=colsNotInGroup(1);
        groupList(i)=currentGroup;
        groupEqns(:)=0;
        groupEqns(a(:,i)~=0)=1;
        for j=2:length(colsNotInGroup)
          colJ=colsNotInGroup(j);
          iaj = find(a(:,colJ));
          if(all(groupEqns(iaj)==0))
            groupList(colJ)=currentGroup;
            groupEqns(iaj) = 1;
          end
        end
        currentGroup = currentGroup + 1;
      end
    end