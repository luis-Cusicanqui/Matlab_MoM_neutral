clear all
close all
%% Define Symbols
M = 8;
f_tilde = sym('f_%d',[M,1],{'real'});
f_tilde(1) = sym('rho',{'positive'});

syms theta u 
x = [f_tilde(1);u;theta;f_tilde(3:end)];
assume(theta,'real');
assume(theta,'positive');
assume(u,'real');


%% Equations
N = M+1;
A = sym(zeros(M+1));

for i=2:N-1
    A(i,i) = x(2);
    A(i,i-1) = x(3);
    A(i,i+1) = i;
end
A(1,1) = u;
A(1,2) = 1;
A(N,N) = u;
A(N,N-1) = x(3);


