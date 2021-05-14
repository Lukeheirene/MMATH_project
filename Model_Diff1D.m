%---------------------------------------------------------------------------
% New 1D DIFFUSION MODEL
%---------------------------------------------------------------------------

m = 0;
x = linspace(-25,25,1000);
t = linspace(0,50,500);

mu_lm; % call parameter tuning function

sol = pdepe(m,@modelpde,@modelic,@modelbc,x,t);
u1  = sol(:,:,1);
u2  = sol(:,:,2);
u3  = sol(:,:,3);

set(0, 'DefaultFigureRenderer', 'opengl');

figure;
surf(x,t,u1,'EdgeColor','none');
title('$c$','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$t$','interpreter','latex');







% --------------------------------------------------------------------------
% NOTE: mu_lm() sets mu as global variable, then acquires
%       cSt, hSt, oSt and sets them as global variables      
% --------------------------------------------------------------------------

function P = mu_lm()

global p;
global lm;

global cSt;
global hSt;
global oSt;

p  = 0.169;
lm  = 0;


R   = Model_SteadyStateSolver(p,lm)

cSt = R(1);
hSt = R(2);
oSt = R(3);

end

function [c, f, s] = modelpde(x,t,u,dudx)

K1 = 4.88566;
K2 = 0.1;
K3 = 0.05;
Ve = 1e-6;
kinfty = 52e-6;
g = 0.5;


% declaring global parameters within function
%--------------------------------------------
global p; 
global lm;

%--------------------------------------------

D  = 0.1;


F  = K1*u(2)*((u(1)^2)/((K2^2) + u(1)^2)) - ((u(1)^2)/(K3^2 + u(1)^2));
G  = ((1/( 1 + (((Ve)/(g*(kinfty*(((p)^4)/(1 + (p)^4)))))*u(1))^4)) - u(2));
O  = u(3);

c = [1; 1; 1];
f = [D; 0; 0] .* dudx;
s = [F; G; O];

end

function u0 = modelic(x)

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

function [pl,ql,pr,qr] = modelbc(xl,ul,xr,ur,t)

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