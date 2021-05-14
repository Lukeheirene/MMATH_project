clear 
clc
close all

%call function ode45 to solve Atri model
[t,CH]=ode45(@Atri,[0 50], [1 1]);
c = CH(:,1);
h = CH(:,2);

%plot solutions
plot(t,c);
hold on;
xlabel('$t$','interpreter','latex');
plot(t,h,'-');
l=legend('$c$ : Ca$^{2+}$ concentration','$h$ : Fraction of active IP$_3$Rs');
set(l, 'interpreter', 'latex')
hold off;

%function setting up Atri model
function A = Atri(t,ch)

%parameters
kf = 16.2 ;
th = 2;
k1 = 0.7;
y = 2;
ky = 0.1;
k2 = 0.7;
b = 0.111;
mu = 0.55;
K1 = (kf*th)/k1;
T = (y*th)/k1;
K = ky/k1;
K2 = k2/k1;

%ODEs
Fch = ((b+ch(1))*mu*K1*ch(2)/(1+ch(1)))-(T*ch(1)/(K+ch(1)));
Gch = (K2^2/(K2^2+ch(1)^2))-ch(2);
A = [Fch; Gch];

end


