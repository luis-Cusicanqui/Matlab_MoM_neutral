clear all
syms x mu sigma
assume(mu,'real')
assume(sigma,'real')
assume(sigma,'positive')
assume(x,'real')
herm = hermite(x,10);
moms = moments(mu,sigma,10);