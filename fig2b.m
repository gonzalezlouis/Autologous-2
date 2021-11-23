clear all

a = 10^-5; %m
eps1 = 0.07;
eps3 = 3*eps1;
w = 2;
alpha = 0.75;
alpha0 = 0;
width = 3*10^-3; %m
rho = logspace(7,12,100);
number = floor(linspace(1,600,20)');
y = floor(logspace(log10(10),log10(600),20)');

%G = load('rhovsanisotropy.csv');
H = readtable('celldensitysimulationcomsol.xlsx');

%cellnumber = G(:,1);
%anisotropy_absorbing = G(:,2); 
%anisotropy_binding = G(:,6);
%anisotropy_alpha0 = G(:,7);
NumberOfCells = H(:,1);
Avg = H(:,7);
Error = H(:,8);
Avg2 = H(:,15);
Error2 = H(:,16);
Avg3 = H(:,23);
Error3 = H(:,24);
Avg4 = H(:,31);
Error4 = H(:,32);

comsol_number = table2array(NumberOfCells);
comsol_anisotropy = table2array(Avg);
comsol_error = table2array(Error);
comsol_anisotropyalpha001 = table2array(Avg2);
comsol_erroralpha001 = table2array(Error2);
comsol_anisotropyalpha01 = table2array(Avg3);
comsol_erroralpha01 = table2array(Error3);
comsol_anisotropyalpha0001 = table2array(Avg4);
comsol_erroralpha0001 = table2array(Error4);

comsol_density = (comsol_number ./ (3000*2000*100)) .* 10^18;

A_1 = (eps1*w/(8*(2 + alpha))) .* (eps1*a)./ ... 
    (eps1.*a + 4*pi*width*(1 + alpha).*rho.*(a^3));

A_3 = (eps3*w/(8*(2 + alpha))) .* (eps3*a)./ ... 
    (eps3.*a + 4*pi*width*(1 + alpha).*rho.*(a^3));

A_alpha0 = (eps1*w/(8*(2 + alpha0))) .* (eps1*a)./ ... 
    (eps1.*a + 4*pi*width*(1 + alpha0).*rho.*(a^3));

A_3alpha0 = (eps3*w/(8*(2 + alpha0))) .* (eps3*a)./ ... 
    (eps3.*a + 4*pi*width*(1 + alpha0).*rho.*(a^3));


%celldensity = (cellnumber ./ (3000*2000*100)) .* 10^18; 
% covnerts the cell density from um^-3 to m^-3

numberdensity = (number ./ (3000*2000*100)) .* 10^18;

figure(1); clf; hold on
lw = 2; lw2 = 1.5; ms = 9; fs = 24; fs2 = 20;

plot([50 50],[1e-4,1e-1],'--','color',[0 .85 0],'linewidth',lw)
plot([250 250],[1e-4,1e-1],'--','color',[0 .45 0],'linewidth',lw)
plot(rho/1e9,A_3,'r','linewidth',lw)
plot(rho/1e9, A_3alpha0, 'k','linewidth',lw);
% plot(celldensity,anisotropy_binding, 'o')
% hold on
errorbar(comsol_density/1e9,comsol_anisotropy, comsol_error, 'ro',...
    'linewidth',lw2,'markersize',ms)
% hold on
% errorbar(comsol_density,comsol_anisotropyalpha001, comsol_erroralpha001, 'o')
% hold on
% errorbar(comsol_density,comsol_anisotropyalpha01, comsol_erroralpha01, '^')
% hold on
errorbar(comsol_density/1e9,comsol_anisotropyalpha0001, comsol_erroralpha0001, 'k^',...
    'linewidth',lw2,'markersize',ms)
set(gca, 'XScale', 'log', 'YScale', 'log','fontsize',fs2);
% ylim([0 0.01])
xlabel('Cell density, $\rho$ (cells/mm$^3$)','fontsize',fs,'interpreter','latex')
ylabel('Anisotropy, $A$','fontsize',fs,'interpreter','latex')
xlim([.5 1000])
ylim([6e-4 5e-2])
%set(gca, 'YDir', 'reverse')
%xline(5e10, '--','low seeding','fontsize',15);
%xline(25e10, '--','high seeding','fontsize',15);
%legend('analytic \alpha = 0.75','comsol \alpha = 0.75',...
%    'analytic \alpha = 0','comsol \alpha = 0.0001',location = 'best')
box on

print(gcf,'-depsc','fig2b.eps')