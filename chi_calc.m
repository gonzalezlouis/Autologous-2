function chi = chi_calc(Y,Z,Y0,Z0,epsilon,kappa,alpha_,beta_)

rho = sqrt((Y-Y0).^2+(Z-Z0).^2);
theta = mod(atan2(Y-Y0,Z-Z0),2*pi); % makes theta go 0 to 2*pi
gamma_ = beta_/(1+alpha_);
w = 1 + 1/kappa - exp(1/kappa)*expint(1/kappa)/kappa^2;
zeta = 1 + 3*kappa + 3*kappa^2;
f = 4 - 4*(2*kappa+1)./rho.^2 + 2*zeta./rho.^3 ...
    + kappa^2*exp(1/kappa).*((rho.^3/kappa^3 ...
    - rho.^2/kappa^2 + 2*rho/kappa - 6).*exp(-rho/kappa) ...
    - rho.^4.*expint(rho/kappa)/kappa^4)./rho.^3;
chi0 = gamma_./rho;
chi1 = gamma_/2*(alpha_/(1+alpha_)./rho - 1 ...
    + cos(theta)/4.*((1-alpha_)*w/(2+alpha_)./rho.^2 + f));
chi = chi0 + epsilon*chi1;
i = find(rho < 1);
chi(i) = 0;
