function beamFreeVibration
E=200e9;
thick = .2;
width = .1;
I = width*thick^3/12;
EI=E*I;
A=thick*width;
L=10;
F=0;
N=0;
xMidPoint = L/2;
rho=7700; % density of steel, kg/m^3
m = 0;
amp=.1;
numElems=20; % even number of elements gives a node in the middle of the beam
elemLength = L/numElems;
numNodes = numElems + 1;
x = linspace(0,L,numNodes);
t=linspace(0, .75, 100);

% use anonymous functions to pass parameters to functions required by pde1dm
pde = @(x,t,u,DuDx) beampde(x,t,u,DuDx, EI, rho, A, F, N);
bc = @(xl,ul,xr,ur,t) beambc(xl,ul,xr,ur,t);
ic = @(x) beamic(x, L, amp);

sol = pde1dm(m,pde,ic,bc,x,t);

% plot the results

midNode = int32(numElems/2 + 1);
uMidPt = sol(:, midNode, 1);
uAnalVibr = analFreeVibr(L, EI, rho*A, amp, x, t);
figure; plot(t, uMidPt, t, uAnalVibr(:,midNode), 'o'); grid;
xlabel Time; ylabel Displacement; 
title('Displacement at Beam Mid-point as a Function of Time');
legend('pde1dm', 'analytical');
figure; plot(x,sol(end,:,1), x, uAnalVibr(end, :), 'o'); grid;
xlabel x; ylabel Displacement; 
title('Beam Displacement at the Final Time');
legend('pde1dm', 'analytical');

end

function [cr,fr,sr] = beampde(x,t,u,DuDx, EI, rho, A, F, N)
cr = [1 rho*A 0]';
fr = [0 -EI*DuDx(3)-N*DuDx(1) -DuDx(1)]';
sr = [u(2) F u(3)]';
end

function u0 = beamic(x, L, amp)
  % half sin wave initial condition
  s = sin(pi*x/L);
  u0 = [amp*s; 0; -amp*(pi/L)^2*s];
end
	
function [pl,ql,pr,qr] = beambc(xl,ul,xr,ur,t)
  pl = [ul(1) ul(2) ul(3)]';
  ql = [0 0 0]';
  pr = [ur(1) ur(2) ur(3)]';
  qr = [0 0 0]';
end

function u=analFreeVibr(L, EI, rhoA, amp, x, t)
 omega=(pi/L)^2*sqrt(EI/rhoA);
 u=amp*cos(omega*t)'*sin(pi*x/L);
 end
  


