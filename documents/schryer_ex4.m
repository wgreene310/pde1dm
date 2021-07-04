function schryer_ex4
n=15;
x = linspace(-pi,pi,n);
t=linspace(0, 3*pi/2.5, 10);
m = 0;
% We need the PDE solution at the right and left ends
% to define the ODE (constraint equation).
xOde = [pi -pi]';
% Set vectorized mode to improve performance
opts.vectorized='on';
[u,vOde] = pde1dm(m,@pdeFunc,@icFunc,@bcFunc,x,t,@odeFunc,@odeIcFunc,xOde,opts);

% Compute analytical solution for comparison.
ua=uAnal(t,x);

figure; plot(x, u(end,:), x, uAnal(t(end), x), 'o'); grid on;
xlabel('x'); ylabel('u'); title('Solution at Final Time');
legend('Numerical', 'Analytical');

figure; plot(t, u(:,end), t, u(:,1), t, uAnal(t, x(1)), 'o'), grid on;
xlabel('Time'); ylabel('u'); title('Solution at End Points');
legend('Right End', 'Left End', 'Analytical');

figure; plot(t, vOde); grid on;
xlabel('Time'); ylabel('v'); title('ODE Variable as a Function of Time');

end

function [c,f,s] = pdeFunc(x,t,u,DuDx,v,vdot)
nx=length(x);
c = ones(1,nx);
f = DuDx;
cx=cos(x);
st=sin(t);
gxt=cx*cos(t)+cx*st+cx.^3*st^3;
s = -u.^3 + gxt;
end

function u0 = icFunc(x)
u0 = uAnal(0, x);
end

function [pl,ql,pr,qr] = bcFunc(xl,ul,xr,ur,t,v,vdot)
% du/dx at right end must equal
% du/dx at the left end. Use the ODE
% variable to enforce this.
pl = v(1);
ql = 1;
pr = v(1);
qr = 1;
end

function f=odeFunc(t,v,vdot,x,u,DuDx)
% The solution at the right end must equal
% the solution at the left end.
f=u(1)-u(2);
end

function v0=odeIcFunc()
v0=0;
end

function u=uAnal(t,x)
u=sin(t)'*cos(x);
end
