%% Copyright (C) 2016, Francesco Faccio <francesco.faccio@mail.polimi.it>
%% Copyright (C) 2016-2017 William H. Greene <w.h.greene@gmail.com>
%%
%% This file is part of Octave.
%%
%% Octave is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or (at
%% your option) any later version.
%%
%% Octave is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn  {} {[@var{y0_new}, @var{yp0_new}] =} decic (@var{fun}, @var{t0}, @var{y0}, @var{fixed_y0}, @var{yp0}, @var{fixed_yp0})
%% @deftypefnx  {} {[@var{y0_new}, @var{yp0_new}] =} decic (@var{fun}, @var{t0}, @var{y0}, @var{fixed_y0}, @var{yp0}, @var{fixed_yp0}, @var{options})
%% @deftypefnx  {} {[@var{y0_new}, @var{yp0_new}, @var{resnorm}] =} decic (@dots{})
%%
%% Compute consistent initial conditions @var{y0_new} and @var{yp0_new}, given
%% initial guesses @var{y0} and @var{yp0}. Choose a maximum of
%% @code{length(@var{y0})} components between @var{fixed_y0} and
%% @var{fixed_yp0} as fixed values.
%%
%% @var{fun} is a function handle. The function must accept three inputs where
%% the first is time @var{t}, the second is a column vector of unknowns @var{y}
%% and the third is a column vector of unknowns @var{yp}.
%%
%% @var{t0} is the initial time such that @code{@var{fun}(@var{t0},
%% @var{y0_new}, @var{yp0_new}) = 0}, specified as a scalar.
%%
%% @var{y0} is a vector used as initial guess for @var{y}.
%%
%% @var{fixed_y0} is a vector which specifies the components of @var{y0} to be
%% hold fixed. Choose a maximum of @code{length(@var{y0})} components between
%% @var{fixed_y0} and @var{fixed_yp0} as fixed values.
%% Set @var{fixed_y0}(i) component to 1 if you want to fix the value of
%% @var{y0}(i).
%% Set @var{fixed_y0}(i) component to 0 if you want to allow the value of
%% @var{y0}(i) to change.
%%
%% @var{yp0} is a vector used as initial guess for @var{yp}.
%%
%% @var{fixed_yp0} is a vector which specifies the components of @var{yp0} to
%% be hold fixed. Choose a maximum of @code{length(@var{yp0})} components
%% between @var{fixed_y0} and @var{fixed_yp0} as fixed values.
%% Set @var{fixed_yp0}(i) component to 1 if you want to fix the value of
%% @var{yp0}(i).
%% Set @var{fixed_yp0}(i) component to 0 if you want to allow the value of
%% @var{yp0}(i) to change.
%%
%% The optional seventh argument @var{options} is a structure array.
%% Use @code{odeset} to generate this structure. The relevant options are
%% @code{RelTol} and @code{AbsTol} which specify the error thresholds used to
%% compute the initial conditions.
%%
%% The function typically returns two outputs. Variable @var{y0_new} is a
%% column vector and contains the consistent initial value of y.  The
%% output @var{yp0_new} is a column vector and contains the consistent initial
%% value of yp.
%%
%% The optional third output @var{resnorm} is the vector of norm of the
%% residuals. If @var{resnorm} is small, @code{decic} has successfully computed
%% the initial conditions. If the value of  @var{resnorm} is large, use
%% @code{RelTol} and @code{AbsTol} to adjust it.
%%
%% Example: Compute initial conditions of @nospell{Robetson}'s equations:
%%
%% @example
%% @group
%% function res = robertsidae(@var{t}, @var{y}, @var{yp})
%% res = [-(@var{yp}(1) + 0.04*@var{y}(1) - 1e4*@var{y}(2)*@var{y}(3));
%%        -(@var{yp}(2) - 0.04*@var{y}(1) + 1e4*@var{y}(2)*@var{y}(3) +
%%        3e7*@var{y}(2)^2);
%%        @var{y}(1) + @var{y}(2) + @var{y}(3) - 1];
%% endfunction
%% @end group
%% [@var{y0_new},@var{yp0_new}] = decic (@@robertsidae, 0, [1; 0; 0], [1; 1; 0],
%% [-1e-4; 1; 0], [0; 0; 0]);
%% @end example
%% @seealso{odeset}
%% @end deftypefn

function [y0_new, yp0_new, resnrm] = decicShampine (odefun, t0, y0, fixed_y0, yp0, ...
  fixed_yp0, options)

if (nargin < 6 || nargin > 7 || nargout > 3)
  print_usage ();
end

%Check input
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


%!function res = rob (t, y, yp)
%! res =[-(yp(1) + 0.04*y(1) - 1e4*y(2)*y(3)); -(yp(2) - 0.04*y(1) +
%!      + 1e4*y(2)*y(3) + 3e7*y(2)^2); y(1) + y(2) + y(3) - 1];
%!endfunction

%!test  % Without options
%! ref1 = [1;0;0];
%! ref2 = [-4e-2; 4e-2; 0];
%! [ynew,ypnew] = decicShampine (@rob,0,[1;0;0],[1;1;0],[23;110;0],[0;0;1]);
%! assert ([ynew(1:end), ypnew(1:end)], [ref1(1:end), ref2(1:end)], 1e-10);
%!test  % With options
%! ref1 = [1;0;0];
%! ref2 = [-4e-2; 4e-2; 0];
%! opt = odeset ("AbsTol", 1e-8, "RelTol", 1e-4);
%! [ynew,ypnew] = decicShampine (@rob,0,[1;0;0],[1;1;0],[23;110;0],[0;0;1],opt);
%! assert ([ynew(1:end), ypnew(1:end)], [ref1(1:end), ref2(1:end)], 1e-5);

%% Test input validation
%!error decicShampine ()
%!error decicShampine (1)
%!error decicShampine (1,2)
%!error decicShampine (1,2,3,4,5,6,7,8)
%!error <FUN must be a valid function handle>
%!  decicShampine (1, 0, [1; 0; 0], [1; 1; 0], [-1e-4; 1; 0], [0; 0; 0]);
%!error <INIT must be a numeric scalar value>
%!  decicShampine (@rob, [1, 1], [1; 0; 0], [1; 1; 0], [-1e-4; 1; 0], [0; 0; 0]);
%!error <length of y0, fixed_y0, yp0 and fixed_yp0 must be equal>
%!  decicShampine (@rob, 0, [0; 0], [1; 1; 0], [-1e-4; 1; 0], [0; 0; 0]);
%!error <y0, fixed_y0, yp0 and fixed_yp0 must be numeric vectors>
%!  decicShampine (@rob, 0, [1; 0; 0], [1; 0],"", [0; 0; 0]);
%!error <y0, fixed_y0, yp0 and fixed_yp0 must be numeric vectors>
%!  decicShampine (@rob, 0, [1; 0; 0], [1; 1; 0], [-1e-4; 1; 0], [0; 0; "1"]);
%!error <fixed_y0 and fixed_yp0 must be boolean vectors>
%!  decicShampine (@rob, 0, [1; 0; 0], [5; 5; 0], [-1e-4; 1; 0], [0; 0; 0]);
%!error <fixed_y0 and fixed_yp0 must be boolean vectors>
%!  decicShampine (@rob, 0, [1; 0; 0], [1; 1; 0], [-1e-4; 1; 0], [0; 4; 0]);
%!error %too many components fixed
%!  decicShampine (@rob, 0, [1; 0; 0], [1; 1; 0], [-1e-4; 1; 0], [0; 1; 1]);


