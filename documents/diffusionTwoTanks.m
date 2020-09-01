function diffusionTwoTanks
L=1;
D=2; % thermal conductivity of conduction region
C=1; % rho*c_p of conduction region
x = linspace(0,L,20);
t=linspace(0,1.5, 80);
m = 0;
% the two ODE will be coupled to the PDE at the
% two ends of the diffusion section
xOde = [0 L]';
icF = @(x) icFunc(x,0);
odeIcF = @() odeIcFunc(0);
pdeF = @(x,t,u,DuDx,v,vdot) pdeFunc(x,t,u,DuDx,v,vdot,D,C);
odeF = @odeFunc;
bcF = @bcFunc;
[u,uode] = pde1dM(m,pdeF,icF,bcF,x,t,odeF, odeIcF,xOde);

figure; plot(x, u(end,:)); grid on;
xlabel('x');
title('Solution at Final Time');
%figure; plot(t, uode, t, vAnal(t)); grid on;
figure; plot(t, uode); grid on;
xlabel('Time'); title('Solution to the two ODE as a Function of Time');
legend('Left Tank', 'Right Tank');

figure; plot(t, u(:,1), t, u(:,end)), grid on;
xlabel('Time');
title('Solution at the two End as a Function of Time');
legend('Left', 'Right');

end

function [c,f,s] = pdeFunc(x,t,u,DuDx,v,vdot,D,C)
c = C;
f = D*DuDx;
s= 0;
end

function u0 = icFunc(x, t0)
% initial temperature in the conduction region
u0=0;
end

function [pl,ql,pr,qr] = bcFunc(xl,ul,xr,ur,t,v,vdot)
% temperatures at the two ends are equal to the temperatures
% in the two tanks
pl = ul-v(1);
ql = 0;
pr = ur-v(2);
qr = 0;
end

function F=odeFunc(t,v,vdot,x,u,DuDx,f,~,~)
F=[vdot(1)+f(1) vdot(2)+f(2)];
end

function v0=odeIcFunc(t0)
% initial temperatures in left and right tanks
v0=[3 0];
end



