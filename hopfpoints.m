clear
clc

%define variables 
%everything with a "_1" are used to calculate the left Hopf points

Gm = 5.71;
K_1 = 46.2857;
K = 0.1429;
b = 0.111;
P = linspace(0,1.4613,1000);
c0 = [0.25 1.5];
c0_1 = [0.2 0.25];

for i=1:length(P)
myfunc = @(c) (Gm*c*(c^3 + c^2 + c + 1))/((K_1)*c^2 + K_1*(K + b)*c + K*K_1*b)- ((((1 + c)^2)*(1 + c^2))/(K_1*(1 - b)))*((Gm/(K + c))*(1 - ((c)/(K + c))) + 1 + P(i)); 
c(i) = fzero(myfunc,c0);
c_1(i) = fzero(myfunc,c0_1);
g(i) = ((((1 + c(i))^2)*(1 + c(i)^2))/(K_1*(1 - b)))*((Gm/(K + c(i)))*(1 - ((c(i))/(K + c(i)))) + 1 + P(i));
g_1(i) = ((((1 + c_1(i))^2)*(1 + c_1(i)^2))/(K_1*(1 - b)))*((Gm/(K + c_1(i)))*(1 - ((c_1(i))/(K + c_1(i)))) + 1 + P(i));
Pp(i) = (Gm*c(i)*(c(i)^3 + c(i)^2 + c(i) + 1)*K_1*(1 - b))/(((K_1)*(c(i)^2) + K_1*(K +b)*c(i) +K_1*K*b)*(1 + c(i)^2)*((1 + c(i))^2)) - (Gm*K)/((K + c(i))^2) - 1;
Pp_1(i) = (Gm*c_1(i)*(c_1(i)^3 + c_1(i)^2 + c_1(i) + 1)*K_1*(1 - b))/(((K_1)*(c_1(i)^2) + K_1*(K +b)*c_1(i) +K_1*K*b)*(1 + c_1(i)^2)*((1 + c_1(i))^2)) - (Gm*K)/((K + c_1(i))^2) - 1;
c0 = c(i);
c0_1 = c_1(i);
end

%variables%
c = linspace(0,5,1000);
mu = ones(1,length(c));
P = ones(1,length(c));

for i=1:length(c)
mu(i) = (Gm*c(i)*(c(i)^3 + c(i)^2 + c(i) + 1))/(K_1*c(i)^2 + K_1*(K + b)*c(i) +K*K_1*b);
 P(i) = (Gm*c(i)*(c(i)^3 + c(i)^2 + c(i) + 1)*K_1*(1 - b))/(((K_1)*(c(i)^2) + K_1*(K +b)*c(i) +K_1*K*b)*(1 + c(i)^2)*((1 + c(i))^2)) - (Gm*K)/((K + c(i))^2) - 1; 
end

%plot%
figure
plot(mu,P,'r o','linewidth',3)
hold on
plot(g,Pp,'b','linewidth',3)
plot(g_1,Pp_1,'b','linewidth',3)
xlabel('\mu')
ylabel('$P$','interpreter','latex')
axis([0.2 0.6 0 1.6])