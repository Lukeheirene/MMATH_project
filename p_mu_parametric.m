clear
clc


%define variables
Gm = 5.71;
K_1 = 46.2857;
K = 0.1429;
b = 0.111;
c = linspace(0,5,1000);
mu = ones(1,length(c));
P = ones(1,length(c));

for i=1:length(c)
mu(i) = (Gm*c(i)*(c(i)^3 + c(i)^2 + c(i) + 1))/(K_1*c(i)^2 + K_1*(K + b)*c(i) +K*K_1*b);
 P(i) = (Gm*c(i)*(c(i)^3 + c(i)^2 + c(i) + 1)*K_1*(1 - b))/(((K_1)*(c(i)^2) + K_1*(K +b)*c(i) +K_1*K*b)*(1 + c(i)^2)*((1 + c(i))^2)) - (Gm*K)/((K + c(i))^2) - 1; 
 %P_1(i) = (-Gm*K)/((K + c(i))^2) + ((mu(i)*K_1)/((1 + c(i)^2)*(1 + c(i))))*((1 - b)/(1 + c(i)) -(2*(b + c(i))*c(i))/(1 + c(i)^2));
end

%ploting of Hopf curve
figure
plot(mu,P)
xlabel('$\mu$','interpreter','latex')
ylabel('$P$','interpreter','latex')
axis([0.25 0.55 0 1.6])
%%
%Atri steady state curve
c_1=linspace(0,0.6,1000);
mu_1 = (Gm.*c_1.*(c_1.^3 + c_1.^2 + c_1 + 1))./(K_1.*c_1.^2 + K_1.*(K + b).*c_1 +K*K_1*b);
plot(mu_1,c_1)
xline(0.2891,'r')
xline(0.28795,'r')
axis([0.27 0.3 0.1 0.6])
xlabel('$\mu$','interpreter','latex')
ylabel('$c$','interpreter','latex')
