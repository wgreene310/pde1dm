function meltingIce
% Melting of a semi-infinite block of ice
n=21;
% Define mesh in parametric space ranging from zero to one.
xi = linspace(0,1,n);

tfinal=3600; % simulate block melting for one hour
nt=30;
t = linspace(0,tfinal,nt);

% physical properties of water
rho=1; % density, gm/cm^3
cp=4.1868; % heat capacity
k=.564e-2; % thermal conductivity
latentHeat=333.4; %  J/g
alpha=k/(rho*cp);
xOde = 1;
TL=25;
Tm=0;
St=cp*(TL-Tm)/latentHeat; % Stefan number

m=0;
opts.vectorized='on'; % use vectorized option to improve performance
% tighten the accuracy tolerances for ODE solver
opts.reltol=1e-5;
opts.abstol=1e-7;
pdeFunc = @(x,t,u,DuDx,X,Xdot) heatpde(x,t,u,DuDx,X,Xdot,alpha);
icFunc = @(x) heatic(x);
bcFunc = @(xl,ul,xr,ur,t) heatbcDir(xl,ul,xr,ur,t,TL,Tm);
odeIcF = @() odeIcFunc();
odef=@(t,X,Xdot,x,u,DuDxi) odeFunc(t,X,Xdot,x,u,DuDxi,latentHeat, rho, k);
[u,X] = pde1dm(m, pdeFunc,icFunc, bcFunc,xi,t,odef, odeIcF,xOde,opts);

[XAnal,uAnal]=stefanOnePhaseAnalSoln(TL,Tm,St,alpha,t,xi);
x=xi'*X';
xAnal=xi'*XAnal';
figure; plot(x(:,end), u(end,:), xAnal(:,end), uAnal(end,:), 'o'); grid; 
title(sprintf('Temperature in the slab at %d seconds', tfinal));
xlabel 'x, cm'; ylabel 'Temperature, degrees-C';
legend('Numerical', 'Analytical');

figure; plot(t, X, t, XAnal, 'o'); ylabel('X, cm'); grid;
xlabel 'Time, seconds'; 
title 'Melting Front Location as a Function of Time';
legend('Numerical', 'Analytical', 'Location','northwest');

end

function [c,f,s] = heatpde(xi,~,u,DuDxi,X,Xdot,alpha)
% x = xi*X(t)
dxDXi=X;
dxDt=Xdot*xi;
% evaluate the PDE at all xi locations (vectorized mode)
nx=length(xi);
c = dxDXi*ones(1,nx);
f = alpha*DuDxi/dxDXi;
s=DuDxi.*dxDt;
end

function u0 = heatic(x)
u0 = 0;
end

function [pl,ql,pr,qr] = heatbcDir(xl,ul,xr,ur,t,TL,Tm)
pl = ul-TL;
ql = 0;
pr = ur-Tm;
qr = 0;
end

function f=odeFunc(t,X,Xdot,x,u,DuDxi,latentHeat, rho, k)
dxDXi=X(1);
f1=rho*latentHeat*Xdot(1)+k*DuDxi/dxDXi;
f=f1;
end

function X0=odeIcFunc()
X0=1e-5; % start with a near-zero initial front position
end

function [frontLocation,temperature]=stefanOnePhaseAnalSoln(TL,Tm,St,alpha,t,xi)
lambda=sqrt(St/2);
f= @(lam) lam*exp(lam^2)*erf(lam)-St/sqrt(pi);
lambda=fzero(f, lambda);
frontLocation=2*lambda*sqrt(alpha*t');
nx=length(xi);
nt=length(t);
t=t(:)';
temperature=zeros(nt,nx);
for i=1:nt
  Xi=frontLocation(i)*xi(:)';
  temperature(i,:)=TL - (TL-Tm)*erf(Xi/(2*sqrt(alpha*t(i))))/erf(lambda);
  temperature(temperature<0) = Tm;
end
end