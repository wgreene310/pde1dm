function heatConduction
% transient heat conduction in a rod
n=12; % number of nodes in x
L=1; % length of the bar in m
x = linspace(0,L,n);
t = linspace(0,2000,30); % number of time points for output of results
m=0; % rectangular Cartesian coordinate system
u = pde1dm(m, @heatpde,@heatic,@heatbc,x,t);

figure; plot(t, u(:,end)); grid on;
xlabel('Time'); ylabel('Temperature');
title('Temperature at right end as a function of time');

figure; plot(x, u(end,:)); grid on;
xlabel('x'); ylabel('Temperature');
title('Temperature along the length at final time');

end

function [c,f,s] = heatpde(x,t,u,DuDx)
rho=8940; % material density
cp=390; % material specific heat
k=385; % material thermal conductivity 
c = rho*cp;
f = k*DuDx;
s = 0;
end

function [pl,ql,pr,qr] = heatbc(xl,Tl,xr,Tr,t)
T0=100; % temperature at left end, degrees C
pl = Tl-T0;
ql = 0;
pr = 0;
qr = 1;
end

function u0 = heatic(x)
T0=100; % temperature at left end, degrees C
if x==0
  u0=T0;
else
  u0=0;
end
end



