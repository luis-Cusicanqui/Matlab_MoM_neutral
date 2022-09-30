clear all
close all
%% Define Symbols
N = 4;
f_tilde = sym('f_%d',[N,1],{'real'});
f_tilde(1) = sym('rho',{'positive'});
f_tilde(2) = 0;
f_tilde(3) = 0;
syms theta u 
assume(theta,'real');
assume(theta,'positive');

%% Equations
B = sym('B',N,'real');
A = sym('A',N);
B = sym(zeros(N));
A = sym(zeros(N));
for i=1:N
    B(i,i) = 1;
    if(i>=2)
        B(i,2) = sqrt(i-1)*f_tilde(i-1);
    end
    if(i>=3)
        B(i,3) = sqrt((i-1)*(i-2))/2*f_tilde(i-2);
    end
    
end


A(1,1) = u;
A(1,2) = f_tilde(1);

A(2,1) = theta;
A(2,2) = u*f_tilde(1);
A(2,3) = f_tilde(1);

A(3,4) = sqrt(3);
A(3,2) = theta * sqrt(2)*f_tilde(1);
A(3,3) = u/2*sqrt(2)*f_tilde(1);
for i=4:N
    A(i,i) = u;
    A(i,i-1) = theta * sqrt(i-1) * ~isequal(f_tilde(i-1),sym(0));
    if(i<N)
    A(i,i+1) = sqrt(i) * ~isequal(f_tilde(i+1),sym(0));
    end
    
    A(i,2) = u*sqrt(i-1)*f_tilde(i-1) + theta*sqrt((i-1)*(i-2))*f_tilde(i-2) ...
        +(i)*f_tilde(i);
    A(i,3) = u*sqrt((i-1)*(i-2))/2*f_tilde(i-2) + sqrt((i-1)*(i-2)*(i-3))/2*theta*f_tilde(i-3) ...
        + sqrt(i-1)*(i)/2*f_tilde(i-1);
end
