function nagD03phfExample
n=10;
L=1;
x = linspace(0,L,n);
t0=1e-4;
t=[t0 .1*2.^[1:5]];
t=linspace(t0, .6, 10);
m = 0;
xOde = L;
icF = @(x) icFunc(x,t0);
odeIcF = @() odeIcFunc(t0);
opts.vectorized='on';
[u,uode] = pde1dM(m, @pdeFunc,icF,@bcFunc,x,t,@odeFunc, odeIcF,xOde,opts);
va=vAnal(t);
ua=uAnal(t,x);

figure; plot(t, uode(:), t, va, 'o'); 
xlabel('Time'); ylabel('v');
legend('Numerical', 'Analytical');
title 'ODE Solution as a Function of Time';
figure; plot(x, u(end,:), x, uAnal(t(end), x), 'o');
xlabel('x'); ylabel('u');
legend('Numerical', 'Analytical');
title 'PDE Solution At the Final Time';
fprintf('Maximum error in ODE valriable=%10.2e\n', max(abs(uode(:)-va(:))));
fprintf('Maximum error in PDE valriable=%10.2e\n', max(abs(u(:)-ua(:))));

end

function [c,f,s] = pdeFunc(x,t,u,DuDx,v,vdot)
nx = length(x);
c = repmat(v(1)^2,1,nx);
f = DuDx;
s = x.*DuDx*v*vdot;
end

function u0 = icFunc(x, t0)
u0 = uAnal(t0, x);
%prtVec('u0', u0);
end

function [pl,ql,pr,qr] = bcFunc(xl,ul,xr,ur,t,v,vdot)
pl = v(1)*exp(t);
ql = 1;
pr = v(1)*vdot(1);
qr = 1;
end

function f=odeFunc(t,v,vdot,x,u,DuDx)
%fprintf('x=%f\n', x);
f=v*u + DuDx + 1 + t - vdot(1);
end

function v0=odeIcFunc(t0)
v0=vAnal(t0);
end

function v=vAnal(t)
  % analytical soluton for ode variable
  v=t;
end

function u=uAnal(t,x)
  % analytical soluton for pde variable
  u=exp(t'*(1-x))-1;
end
