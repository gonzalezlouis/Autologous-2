clear all

muaD = pi*(sqrt(17)-1);
nua2D = 2/3;
epsilon = 1/5;
f = 1/(1+muaD/8/pi);

N = 1e3;
x = linspace(0,.3,N); % c_0*a^3
A1 = epsilon*nua2D/8./(nua2D + 4*pi*x);
A2 = f*epsilon*(nua2D - muaD*x)/8./(nua2D + 4*pi*x);

figure(1); clf
lw = 2; lw2 = 1;
fs = 24; fs2 = 20;
h = plot([min(x) max(x)],[0 0],'k-',x,A1,'k-',x,A2,'r-');
set(h(2:3),'linewidth',lw)
set(gca,'fontsize',fs2,'ytick',-.01:.01:.02)
xlabel('Background concentration, $c_0a^3$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Anisotropy, $A$','fontsize',fs,'interpreter','latex')
legend(h(2:3),{'Binding','Absorption'},'location','sw')
print(gcf,'-depsc','fig4.eps')