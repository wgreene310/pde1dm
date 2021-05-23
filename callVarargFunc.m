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

function varargout=callVarargFunc(func, args)
numExpectedArgs=nargin(func);
if numExpectedArgs>length(args)
  fn = func2str(func);
  error(['Function \"' fn '\" is expecting %d arguments but is being called with' ...
    ' only %d arguments.'], numExpectedArgs, length(args));
end
[varargout{1:nargout}]=func(args{1:numExpectedArgs});
end