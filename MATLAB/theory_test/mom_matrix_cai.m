% clear all
close all
M = 4;
f_tilde = sym('f_%d',[M,1],{'real'});
f_tilde(1) = sym('rho',{'positive'});

syms theta u 
x = [f_tilde(1);u;theta;f_tilde(3:end)];
assume(theta,'real');
assume(theta,'positive');
assume(u,'real')
N= M+1;
%% Equations

B = sym(zeros(M+1));
A_orig = sym(zeros(M+1));
A = sym(zeros(M+1));

for i=1:M+1
    B(i,i) = 1;
    if(i>=2)
        B(i,2) = x(i-1)*~isequal((i-3)*(i-4),0);
    end
    if(i>=3)
        B(i,3) = 1/2*x(i-2)*~isequal((i-4)*(i-5),0);
    end
    
end

% 
% A(1,1) = x(2);
% A(1,2) = x(1);
% 
% A(2,1) = theta;
% A(2,2) = x(2)*x(1);
% A(2,3) = x(1);
% 
% A(3,4) = 3;
% A(3,2) = x(3) *x(1);
% A(3,3) = x(2)/2*x(1);
% for i=4:M+1
%     A(i,i) = x(2);
%     A(i,i-1) = x(3) * ~isequal((i-4),0);
%     if(i<M+1)
%     A(i,i+1) = i;
%     end
%     
%     A(i,2) = x(2)*x(i-1)*~isequal((i-4),0) + x(3)*x(i-2)*~isequal((i-4)*(i-5),0)...
%         +(i)*x(i);
%     A(i,3) = x(2)/2*x(i-2)*~isequal((i-4)*(i-5),0) ...
%         + 1/2*x(3)*x(i-3)*~isequal((i-5)*(i-6),0) ...
%         +(i)/2*x(i-1)*~isequal((i-4),0);
% end
% 
% S = inv(B);
% A = S*A;

A_orig(1,1) = x(2);
A_orig(1,2) = x(1);
A_orig(2,1) = x(3)/x(1);
A_orig(2,2) = x(2);
A_orig(2,3) = 1;

A_orig(3,2) = A_orig(3,2)+2*x(3);
A_orig(3,3) = A_orig(3,3)+ x(2);
A_orig(3,4) = A_orig(3,4)+6/x(1);
for i=4:N
    A_orig(i,i) = x(2);
    A_orig(i,i-1) = x(3) * ~isequal((i-4),0);
    if(i<N)
        A_orig(i,i+1) = i;
    end
    A_orig(i,1) = -x(3)/x(1)*x(i-1)*~isequal((i-4),0);
    A_orig(i,2) = (i)*x(i);
    A_orig(i,3) = A_orig(i,3)+1/2*x(3)*x(i-3)*~isequal((i-5)*(i-6),0) ...
        +(i-2)/2*x(i-1)*~isequal((i-4),0);
    A_orig(i,4) =A_orig(i,4) -3/x(1)*x(i-2)*~isequal((i-4)*(i-5),0);
end


A_orig(M+1,2) = 0;
A_orig(M+1,3) = -x(M) + x(M-2)*x(3)/2;


N=M+1;
A(1,1) = x(2);
A(1,2) = x(1);
A(2,1) = x(3)/x(1);
A(2,2) = x(2);
A(2,3) = 1;

A(3,2) = A(3,2)+2*x(3);
A(3,3) = A(3,3)+ x(2);
A(3,4) = A(3,4)+6/x(1);
for i=4:N
    A(i,i) = x(2);
    A(i,i-1) = x(3) * ~isequal((i-4),0);
    if(i<N)
        A(i,i+1) = i;
    end
    A(i,1) = -x(3)/x(1)*x(i-1)*~isequal((i-4),0);
    A(i,2) = (i)*x(i);
    A(i,3) = A(i,3)+1/2*x(3)*x(i-3)*~isequal((i-5)*(i-6),0) ...
        +(i-2)/2*x(i-1)*~isequal((i-4),0);
    A(i,4) =A(i,4) -3/x(1)*x(i-2)*~isequal((i-4)*(i-5),0);
end


A(N,2) = 0;
A(N,3) = -x(N-1) + x(N-3)*~isequal((N-5)*(N-6),0)*x(3)/2;
