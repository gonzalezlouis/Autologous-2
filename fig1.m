clear all
fd = '../fig/';

% Y = y/a;
% Z = z/a;
% rho = r/a;
% u = v/v_0;

% Analytic parameters
a = 10; % um
D = 150; % um^2/s
v0 = 15; % um/s
epsilon = v0*a/D;
kappa = .1;
alpha_ = (sqrt(17)-1)/4;
nu = 1; % 1/s
beta = nu/(4*pi*a^2);
beta_ = beta*a^4/D;

% Define the space
Ly = 10/2;
Lz = 10/2;
N = 1e3;
Y = linspace(-Ly,Ly,N);
Z = linspace(-Lz,Lz,N);
Ymat = Y'*ones(1,length(Z));
Zmat = ones(length(Y),1)*Z;

Y0s = [0 6 7 -7 -6]/2;
Z0s = [0 -4 5 -7 6]/2;
chi = 0;
for i = 1:length(Y0s)
    chi_i = chi_calc(Ymat,Zmat,Y0s(i),Z0s(i),epsilon,kappa,alpha_,beta_);
    chi_i(chi_i < 0) = 0;
    chi = chi + chi_i;
end
chi(chi < 0) = 0;

figure(1); clf
subplot(1,10,1:9)
hold on
imagesc(Z,Y,chi)
%contour(Z,Y,chi,10,'linewidth',2)
for i = 1:length(Y0s)
    Z1 = linspace(-1,1,N)+Z0s(i);
    Y1 = sqrt(1-(Z1-Z0s(i)).^2)+Y0s(i);
    Z2 = Z1(end:-1:1);
    Y2 = -sqrt(1-(Z2-Z0s(i)).^2)+Y0s(i);
    fill([Z1 Z2],[Y1 Y2],'w','linewidth',1)
end
xlim([min(Z) max(Z)])
ylim([min(Y) max(Y)])
xlabel('z/a')
ylabel('y/a')
title(['\chi = \chi_0 + \epsilon\chi_1 (\epsilon = ' ...
    num2str(round(100*epsilon)/100) ')'])
colormap(bluewhitered)
colorbar
set(gca,'ydir','normal','layer','top')
box on
print(gcf,'-depsc','fig1.eps')

