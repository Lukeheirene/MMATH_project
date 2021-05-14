clear
clc

%parameter values%

Gm = 5.71;
K_1 = 46.2857;
K = 0.1429;
b = 0.111;
P = [0 0.5 1 1.5 2];

%variables%

c = linspace(0,5,1000);
mu_1 = ones(1,length(c)); %defines mu values of calcium steady state in 1D Atri model
mu_2 = ones(1,length(c)); %defines mu values corresponding to Hopf bifurcations(Tr=0) of Atri diffusion model
F = ones(1,length(c));

 figure
for j=1:length(P)

    for i=1:length(c)

        mu_1(i) = (Gm*c(i)*(c(i)^3 + c(i)^2 + c(i) + 1))/(K_1*c(i)^2 + K_1*(K + b)*c(i) +K*K_1*b);

        mu_2(i) = ((((1 + c(i))^2)*(1 + c(i)^2))/(K_1*(1 - b)))*((Gm/(K + c(i)))*(1 - ((c(i))/(K + c(i)))) + 1 + P(j));
        
    end
    plot(mu_2,c)
    hold on 
    
end

plot(mu_1,c)
xlabel('$\mu$','interpreter','latex')
ylabel('$c$','interpreter','latex')
axis([0 1 0 3])
legend('$P=0$','$P=0.5$','$P=1$','$P=1.5$','$P=2$','Ca$^{2+}$ steady state','interpreter','latex')
