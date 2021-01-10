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

function prtShortVec( mat, title, fid, form, pruneTol)
if(nargin < 3 || isempty(fid))
  fid=1;
end
if(nargin < 4)
  form = '%g ';
end
if nargin > 4
  mat(abs(mat)<pruneTol) = 0;
end

mat = mat(:);
n=length(mat);
fprintf(fid, '%s(%d): ', title, n);
for i=1:n
  fprintf(fid, form, mat(i));
end
fprintf(fid, '\n');
end

