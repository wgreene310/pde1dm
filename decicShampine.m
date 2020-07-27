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
% Copyright (C) 2016-2020 William H. Greene

function [y0_new, yp0_new, resnrm] = decicShampine (odefun, t0, y0, fixed_y0, yp0, ...
  fixed_yp0, options)

if (nargin < 6 || nargin > 7)
  error('decic:nrhs', ...
    'Illegal number of input arguments passed to decic.');
end

p=inputParser;
validHandle = @(x) isa(x, 'function_handle');
if (~ isa (odefun, "function_handle"))
  error ('Octave:invalid-input-arg', ...
    'decic: FUN must be a valid function handle');
end

if (~ isnumeric (t0) || numel (t0) ~= 1)
  error ('Octave:invalid-input-arg', ...
    'decic: INIT must be a numeric scalar value');
end

y0 = y0(:);
yp0=yp0(:);

if (~ isnumeric (y0) || ~ isvector (y0) || ~ isnumeric (fixed_y0) || ...
    ~ isvector (fixed_y0) || ~ isnumeric (yp0) || ~ isvector (yp0)|| ...
    ~ isnumeric (fixed_yp0) || ~ isvector (fixed_yp0))
  error ('Octave:invalid-input-arg', ...
    'decic: y0, fixed_y0, yp0 and fixed_yp0 must be numeric vectors');
end

if (~ isequal (numel (y0), numel (fixed_y0), numel (yp0), ...
    numel (fixed_yp0)))
  error ('Octave:invalid-input-arg', ...
    'decic: length of y0, fixed_y0, yp0 and fixed_yp0 must be equal');
end

for i = 1:numel (y0)
  if (~ (fixed_y0(i) == 0 || fixed_y0(i) == 1) || ~ (fixed_yp0(i) == 0 ...
      || fixed_yp0(i) == 1))
    error ('Octave:invalid-input-arg', ...
      'decic: fixed_y0 and fixed_yp0 must be boolean vectors');
  end
end

n  = numel (y0);
free_y0 = ~fixed_y0;
nl = sum (free_y0);
anyFreeY0 = nl>0;
free_yp0 = ~fixed_yp0;
nu = sum (free_yp0);
anyFreeYp0 = nu>0;

if (n - nl - nu > 0)
  error ('Octave:invalid-input-arg', ...
    'decic: you cannot fix more than length(y0) components');
end

% make sure they are column vectors
y0 = y0(:); yp0 = yp0(:);

%Set default value
relTol = 1e-3;
absTol   = 1e-6;

haveJac = 0;
icdiag=0;
%Check AbsTol, RelTol, Jacobian
if (nargin == 7 && ~ isempty(options))
  if (isfield(options, 'AbsTol'))
    if (~ isscalar (options.AbsTol))
      error ('Octave:invalid-input-arg', ...
        'decic: AbsTol must be a scalar value');
    else
      absTol = options.AbsTol;
    end
  end
  
  if (isfield(options, 'RelTol'))
    if (~ isscalar (options.RelTol))
      error ('Octave:invalid-input-arg', ...
        'decic: RelTol must be a scalar value');
    else
      relTol = options.RelTol;
    end
  end
  
  icdiag=getfieldi(options,'ICDiagnostics');
  if(~icdiag)
    icdiag = 0;
  end
  
  if(isfield(options, 'Jacobian') && ~isempty(options.Jacobian))
    if(iscell(options.Jacobian) && length(options.Jacobian)==2)
      dfDy = options.Jacobian{1};
      dfDyp = options.Jacobian{2};
      checkJacSize(dfDy, dfDyp, n);
      [md,nd]=size(dfDy); [mp,np]=size(dfDy);
      haveJac = 1;
    elseif(isa(options.Jacobian,'function_handle'))
      haveJac = 2;
    else
      error ('Octave:invalid-input-arg', ...
        'decic: Invalid value for Jacobian option');
    end
  end
else
  options = [];
end

%fprintf('anyFreeY0=%d, anyFreeYp0=%d\n', anyFreeY0, anyFreeYp0);

maxIter = 10;
it = 0;
y0_new = y0; yp0_new = yp0;
if(icdiag)
  prtShortVec(y0, 'y0');
  prtShortVec(yp0, 'yp0');
end
while(it <= maxIter)
  res = odefun(t0,y0_new, yp0_new);
  resnrm=norm(res);
  maxRes = max(max(relTol*abs(res),absTol));
  if(icdiag > 1)
    fprintf('iteration=%d, maxres=%12.3e\n', it, maxRes);
  end
  if(maxRes <= absTol)
    return;
  end
  res=res(:);
  
  dfDy = []; dfDyp = [];
  if(haveJac == 0)
    % calculate by forward difference
    if(anyFreeY0)
      dfDy=getDfDy(odefun, t0, y0_new, yp0_new, res, options);
    end
    if(anyFreeYp0)
      dfDyp=getDfDyp(odefun, t0, y0_new, yp0_new, res, options);
    end
  elseif(haveJac == 2)
    [dfDy, dfDyp] = options.Jacobian(t0, y0_new, yp0_new);
    checkJacSize(dfDy, dfDyp, n);
  end
  %dfDy
  %dfDyp
  
  % remove columns for the components of yp0 that are fixed
  dfDyp = dfDyp(:, free_yp0);
  % remove columns for the components of y0 that are fixed
  dfDy = dfDy(:, free_y0);
  
  if (~anyFreeY0)
    % must do QR on full matrix to get full rank-revealing algorithm
    [Q,R,E] = qr(full(dfDyp));
    rr=rankR(R);
    dd=1:rr;
    d = -Q'*res;
    yp0_new =   yp0_new + E(:,dd)*(R(dd,dd)\d(dd));
  elseif (~anyFreeYp0)
    % must do QR on full matrix to get full rank-revealing algorithm
    [Q,R,E] = qr(full(dfDy));
    rr=rankR(R);
    dd=1:rr;
    d = -Q'*res;
    y0_new = y0_new  + E(:,dd)*(R(dd,dd)\d(dd));
  else
    % must do QR on full matrix to get full rank-revealing algorithm
    [Q,R,E] = qr(full(dfDyp));
    algdofs = abs(diag(R))<100*eps;
    %algdofs'
    anyAlgDofs = any(algdofs);
    difdofs=~algdofs;
    R11 = R(difdofs, difdofs);
    %rank(R11)
    d = -Q'*res;
    S = Q'*dfDy;
    if(anyAlgDofs)
      S2122 = S(algdofs,:);
      %S2122
      if(0)
        w=S2122\d(algdofs);
      else
        [Qs,Rs,Es] = qr(S2122);
        % FIXME: need to check rank here
        w = Es*(Rs\(Qs'*d(algdofs)));
      end
      y0_new(free_y0) = y0_new(free_y0) + w;
      w1p = R11\(d(difdofs) - S(difdofs,:)*w);
    else
      w1p = R11\d(difdofs);
    end
    wp = zeros(nu,1);
    wp(difdofs) = w1p;
    yp0_new(free_yp0) = yp0_new(free_yp0) + E*wp;
  end
  it = it + 1;
end

warning ('decic: Failed to obtain a converged set of consistent initial conditions.',...
  ' This might cause the ODE to DAE solver to fail in the first step.');

end

function dfDy=getDfDy(odefun, t, y, yp, r, options)
n = length(y);
dfDy = zeros(n,n);
sqrtEps = sqrt(eps);
for i=1:n
  ysave = y(i);
  delta = sqrtEps*max(abs(y(i)), 1);
  y(i) =   y(i) + delta;
  rp=odefun(t, y, yp);
  dfDy(:,i) = (rp-r)/delta;
  y(i) = ysave;
end
end

function dfDyp=getDfDyp(odefun, t, y, yp, r, options)
n = length(y);
dfDyp = zeros(n,n);
sqrtEps = sqrt(eps);
for i=1:n
  ypsave = yp(i);
  delta = sqrtEps*max(abs(yp(i)), 1);
  yp(i) = yp(i) + delta;
  rp=odefun(t, y, yp);
  dfDyp(:,i) = (rp-r)/delta;
  yp(i) = ypsave;
end
end

function checkJacSize(dfDy, dfDyp, n)
[md,nd]=size(dfDy); [mp,np]=size(dfDy);
if(md~=n || nd~=n || mp~=n || np~=n)
  msg=sprintf('decic: Jacobian matrices must be %d x %d', n, n);
  error ('Octave:invalid-input-arg', msg);
end
end

function rank=rankR(R)
adr=abs(diag(R));
tol=eps*max(adr);
rank=nnz(adr>tol);
end

function value = getfieldi(S,field)
names   = fieldnames(S);
isField = strcmpi(field,names);

if any(isField)
  value = S.(names{isField});
else
  value = [];
end
end
