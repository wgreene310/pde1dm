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

function prtMat( mat, title, fid, format, pruneTol)
if(nargin < 3 || isempty(fid))
  fid=1;
end
if(nargin<4)
  format = '%g,';
end
if nargin > 4
  mat(abs(mat)<pruneTol) = 0;
end
[m,n]=size(mat);
fprintf(fid, '%s(%d,%d)\n', title, m, n);
for r=1:m
  for c=1:n
    fprintf(fid, format, full(mat(r,c)));
  end
  fprintf(fid, '\n');
end
end

