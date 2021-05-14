%---------------------------------------------------------------------------
% ATRI MODEL WITH MECHANICS - 1D DIFFUSION MODEL
%---------------------------------------------------------------------------

m = 0;
x = linspace(-25,25,1000);
t = linspace(0,50,200);

mu_lm; % call parameter tuning function

sol = pdepe(m,@atripde,@atriic,@atribc,x,t);
u1  = sol(:,:,1);
u2  = sol(:,:,2);
u3  = sol(:,:,3);

% Makes plot curve smoother and sets lines visible on legend 
% set(0, 'DefaultFigureRenderer', 'painters');

set(0, 'DefaultFigureRenderer', 'opengl');


figure;
surf(x,t,u1,'EdgeColor','none');
title('$c$','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$t$','interpreter','latex');
supersizeme(1.7)









function P = mu_lm()

global mu;
global lm;
global cSt;
global hSt;
global oSt;

mu  = 0.45;
lm  = 0;


R   = AtriMech_SteadyStateSolver(mu,lm)


cSt = R(1);
hSt = R(2);
oSt = R(3);

end

function [c, f, s] = atripde(x,t,u,dudx)

kf = 16.2;
th = 2;
k1 = 0.7;
y  = 2;
ky = 0.1;
k2 = 0.7;
b  = 0.111;

% declaring global parameters within function
%--------------------------------------------
global mu; 
global lm;
%--------------------------------------------

K1 = (kf*th)/k1;
T  = (y*th)/k1;
K  = ky/k1;
K2 = k2/k1;
D  = 10 ;


F  = ((b+u(1))*mu*K1*u(2)/(1+u(1)))-(T*u(1)/(K+u(1)));
G  = (K2^2/(K2^2+u(1)^2))-u(2);


c = [1; 1; 1];
f = [D; 0; 0] .* dudx;
s = [F; G; 0];

end

function u0 = atriic(x)

% declaring global parameters within function
%--------------------------------------------
global cSt;
global hSt;
global oSt;
%--------------------------------------------

S   = 10; 
sig = 0.01;

u0 = [ cSt+4*exp(-0.5*(x^2)/(S*sig)^2); hSt; oSt];

end

function [pl,ql,pr,qr] = atribc(xl,ul,xr,ur,t)

% declaring global parameters within function
%--------------------------------------------
global cSt;
global hSt;
global oSt;
%--------------------------------------------


% NO FLUX CONDITION
pl = [0; 0; 0];
ql = [1; 1; 1];
pr = [0; 0; 0];
qr = [1; 1; 1];



end