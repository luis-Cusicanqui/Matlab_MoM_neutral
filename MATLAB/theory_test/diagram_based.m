% clear all
close all
%% Define Symbols
M = 4;
f_tilde = sym('f_%d',[M,1],{'real'});
f_tilde(1) = sym('rho',{'positive'});

syms theta u 
x = [f_tilde(1);u;theta;f_tilde(3:end)];
assume(theta,'real');
assume(theta,'positive');
assume(u,'real')

%% Equations

B = sym(zeros(M+1));
A = sym(zeros(M+1));
for i=1:M+1
    B(i,i) = 1;
    if(i>=2)
        B(i,2) = x(i-1)*~isequal((i-3)*(i-4),0);
    end
    if(i>=3)
        B(i,3) = 1/2*x(i-2)*~isequal((i-4)*(i-5),0);
    end
%     B(i,3) = B(i,3)+ x(i)/2*~isequal((i-3)*(i-2),0);
end


A(1,1) = x(2);
A(1,2) = x(1);

A(2,1) = theta;
A(2,2) = x(2)*x(1);
A(2,3) = x(1);

A(3,4) = 3;
A(3,2) = x(3) *x(1);
A(3,3) = x(2)/2*x(1);
for i=4:M+1
    A(i,i) = x(2);
    A(i,i-1) = x(3) * ~isequal((i-4),0);
    if(i<M+1)
    A(i,i+1) = i;
    end
    
    A(i,2) = x(2)*x(i-1)*~isequal((i-4),0) + x(3)*x(i-2)*~isequal((i-4)*(i-5),0)...
        +(i)*x(i);
    A(i,3) = x(2)/2*x(i-2)*~isequal((i-4)*(i-5),0) ...
        + 1/2*x(3)*x(i-3)*~isequal((i-5)*(i-6),0) ...
        +(i)/2*x(i-1)*~isequal((i-4),0);
end
M_grad = A;
B_grad = B;

M_jul = M_grad;
B_jul = B_grad;

M_hme = M_grad;
B_hme = B_grad;


B_jul(end,3) = B_grad(end,3) - (M+1)/2/theta*x(end);

M_jul(end-1,3) = M_jul(end-1,3) - M*(M+1)/2/theta*x(end);
M_jul(end,3) = M_jul(end,3) - (M+1)/2*x(end-1);
M_jul(end,2) = M_jul(end,2) - (M+1)*x(end);
M_jul(end,3) = M_jul(end,3) - u*(M+1)/2/theta * x(end);

M_hme(end,2) = M_hme(end,2) - (M+1)*x(end);
M_hme(end,3) = M_hme(end,3) - (M+1)/2*x(end-1);

A_hme = inv(B_hme)*M_hme;

A_jul = simplify(inv(B_jul)*M_jul);
