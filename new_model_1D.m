clear 
clc
close all

%call function ode45 to solve new model
[t,CN]=ode45(@Model,[0 100],[1 1]);
c = CN(:,1);
h = CN(:,2);

%plot solutions
plot(t,c);
hold on;
xlabel('$t$','interpreter','latex');
plot(t,h,'-');
l=legend('$c$ : Ca$^{2+}$ concentration','$h$ : Fraction of active IP$_3$Rs');
set(l, 'interpreter', 'latex')
hold off;

%function setting up the new model
function M = Model(t,cn)

%parameters
K1 = 4.88566;
K2 = 0.1;
K3 = 0.05;
p=0.3;
kinfty = 52e-6;
Ve = 1e-6;
g=0.5;
kinh = kinfty*(((p)^4)/(1 + (p)^4));
K4 = (Ve)/(g*kinh);

%ODEs
Fcn = K1*cn(2)*((cn(1)^2)/((K2^2) + cn(1)^2)) - ((cn(1)^2)/(K3^2 + cn(1)^2));
Gcn = ((1/( 1 + (K4*cn(1))^4)) - cn(2));

M = [Fcn; Gcn];
end
