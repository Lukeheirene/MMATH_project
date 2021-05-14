clear all
clc
close all

%define variables
pmax=0.81;
Hact=2;
Hinh=3.9;
Kip3=50e-9;
Hip=4;
Kinf=52e-6;
p=linspace(5e-9,100e-9,5000);
p_1=linspace(10e-9,100e-6,5000);
c=linspace(10e-9,100e-6,5000);
Kact = 210e-9;
Kinh_0 = 42e-6;
Kinh_1 = 54e-6;
Kinh_2 = 59e-6;
Kinh_3 = 9.5e-6;
Kinh_4 = 210e-9;
Kinh_5 = 160e-9;

%Plot of biphasic Hill equation
figure
P0_1 = pmax.*((1 + ((Kact)./(c)).^Hact).^(-1)).*((1 + ((c)./Kinh_0).^Hinh).^(-1));
y = semilogx(c,P0_1);
y.LineWidth = 2;
hold on
P0_1 = pmax.*((1 + ((Kact)./(c)).^Hact).^(-1)).*((1 + ((c)./Kinh_1).^Hinh).^(-1));
y = semilogx(c,P0_1);
y.LineWidth = 2;
P0_1 = pmax.*((1 + ((Kact)./(c)).^Hact).^(-1)).*((1 + ((c)./Kinh_2).^Hinh).^(-1));
y = semilogx(c,P0_1);
y.LineWidth = 2;
P0_1 = pmax.*((1 + ((Kact)./(c)).^Hact).^(-1)).*((1 + ((c)./Kinh_3).^Hinh).^(-1));
y = semilogx(c,P0_1);
y.LineWidth = 2;
P0_1 = pmax.*((1 + ((Kact)./(c)).^Hact).^(-1)).*((1 + ((c)./Kinh_4).^Hinh).^(-1));
y = semilogx(c,P0_1);
y.LineWidth = 2;
P0_1 = pmax.*((1 + ((Kact)./(c)).^Hact).^(-1)).*((1 + ((c)./Kinh_5).^Hinh).^(-1));
y = semilogx(c,P0_1);
y.LineWidth = 2;
xlabel('[Ca$^{2+}$]$_i$','interpreter','latex')
ylabel('$P_o$','interpreter','latex')
legend(' $K_{inh}$ = 42$\mu$ M', ' $K_{inh}$ = 54$\mu$ M', ' $K_{inh}$ = 59$\mu$ M',' $K_{inh}$ = 9.5$\mu$ M', ' $K_{inh}$ = 210 nM', ' $K_{inh}$ = 160 nM','interpreter','latex')

%Plot of Hill equation for Kinh

Kinh = Kinf*((1 + (Kip3./p_1).^Hip ).^(-1));
x = loglog(p_1,Kinh);
xlabel('[IP$_3$]','interpreter','latex')
ylabel('$K_{inh}$','interpreter','latex')
x.LineWidth = 2;

%Full dependance of open probability on IP3 and Ca plot

figure
[P,C]= meshgrid(p,c);
Kinh = Kinf*((1 + (Kip3./P).^Hip ).^(-1));
P0 = pmax.*((1 + ((Kact)./(C)).^Hact).^(-1)).*((1 + ((C)./Kinh).^Hinh).^(-1));
surf(C,P,P0)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('[Ca$^{2+}$]$_i$','interpreter','latex')
ylabel('[IP$_3$]','interpreter','latex')
zlabel('$P_o$','interpreter','latex')
shading interp 

