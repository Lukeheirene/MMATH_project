clear
close all
clc

%variables

K1 = 4.88566;
K2 = 0.1;
K3 = 0.05;
Hact = 2;
kinfty = 52e-6;
Ve = 1e-6;
c = linspace(0,0.5,1000);
g = 0.5;

%steady state values

K4 = (((K1-1).*(c.^2) - (K2^2) + (K1*(K3^2)))./((c.^6) + (K2.*(c.^4)))).^(1/4);
p = ((Ve)./(((g*kinfty).*K4) - Ve)).^(1/4);

%steady state curve

plot(p,c)
xlabel('$p$','interpreter','latex')
ylabel('$c$','interpreter','latex')
axis([0 0.5 0 c(end)])

