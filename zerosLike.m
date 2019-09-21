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
%   Copyright (C) 2016-2019 William H. Greene

function z=zerosLike(varargin)
persistent isOctave
if isempty(isOctave)
  isOctave = exist('OCTAVE_VERSION', 'builtin');
end
x=varargin{end};
%z = builtin('repmat', x(1), m, n);
if(isOctave)
  z = x(1)*zeros(varargin{1:end-1});
else
  z = zeros(varargin{1:end-1}, 'like', x);
end
end