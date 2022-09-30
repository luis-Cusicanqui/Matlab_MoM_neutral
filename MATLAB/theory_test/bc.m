clear all
syms x theta rho
assume(rho,'real');
assume(theta,'real');
assume(theta,'positive');

her = hermite(x,10);
f_max = rho/sqrt(2*sym(pi)*theta)*exp(-x^2/(2*theta));
f_mom = 1/sqrt(2*sym(pi))*exp(-x^2/2);

%i=1;
int(x*her(1:5)*f_mom,x,-inf,0)

%i=3
i=3
int(x^i*her(1:5)*f_mom,x,-inf,0)

%i=5
i=5
int(x^i*her(1:5)*f_mom,x,-inf,0)