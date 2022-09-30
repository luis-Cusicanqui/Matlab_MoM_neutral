clear all
close all
%% Define Symbols
M = 5;
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
A_grad = A;
S = inv(B);
P = S*A;
A_real = P;
A_jul(1,1) = x(2);
A_jul(1,2) = x(1);
A_jul(2,1) = x(3)/x(1);
A_jul(2,2) = x(2);
A_jul(2,3) = 1;

A_jul(3,2) = 2*x(3);
A_jul(3,3) = x(2);
A_jul(3,4) = 6/x(1);

N=M+1;
for i=4:N
    A_jul(i,i) = x(2);
    A_jul(i,i-1) = x(3) * ~isequal((i-4),0);
    if(i<N)
    A_jul(i,i+1) = i;
    end
    
    A_jul(i,1) = A_jul(i,1)-x(3)/x(1)*x(i-1)*~isequal((i-4),0);
    A_jul(i,2) =A_jul(i,2)+ (i)*x(i);
    A_jul(i,3) = A_jul(i,3)+ 1/2*x(3)*x(i-3)*~isequal((i-5)*(i-6),0) ...
        +(i-2)/2*x(i-1)*~isequal((i-4),0);
    A_jul(i,4) = A_jul(i,4) -3/x(1)*x(i-2)*~isequal((i-4)*(i-5),0);
end

i=N-1;
A_jul(i,3) = (i-2)/2*(x(i-1)*~isequal((i-3)*(i-4),0))+x(3)*x(i-3)/2*~isequal((i-5)*(i-6),0)...
    - (i+1)*(i)/2*x(i+1)/x(3);
i=N;
A_jul(i,3) = - x(i-1)+ x(i-3)*~isequal((i-5)*(i-6),0)*x(3)/2;
A_jul(i,4) = -3*x(i-2)*~isequal((i-4)*(i-5),0)/x(1) + 3*(i)*x(i)/x(1)/x(3);



% N=length(x);
% A = sym(zeros(N,N));
% for i=4:N
%     A(i,i) = x(2);
%     A(i,i-1) = x(3) * ~isequal((i-4),0);
%     if(i<N)
%     A(i,i+1) = i;
%     end
%     
%     A(i,1) = -x(3)/x(1)*x(i-1)*~isequal((i-4),0);
%     A(i,2) = (i)*x(i);
%     A(i,3) = A(i,3)+1/2*x(3)*x(i-3)*~isequal((i-5)*(i-6),0) ...
%         +(i-2)/2*x(i-1)*~isequal((i-4),0);
%     A(i,4) =A(i,4) -3/x(1)*x(i-2)*~isequal((i-4)*(i-5),0);
% end
% 
% % QMBE
% 
% i=N-1;
% A(i,3) = (i-2)/2*(x(i-1)*~isequal((i-3)*(i-4),0))+x(3)*x(i-3)/2*~isequal((i-5)*(i-6),0)...
%     - (i+1)*(i)/2*x(i+1)/x(3);
% i=N;
% % A(i,2) = 0;
% A(i,3) = - x(i-1)+ x(i-3)*~isequal((i-5)*(i-6),0)*x(3)/2;
% A(i,4) = -3*x(i-2)*~isequal((i-4)*(i-5),0)/x(1) + 3*(i)*x(i)/x(1)/x(3);
% if(M==4)
%     A(i,4) = A(i,4)+ x(3);
% end
% 
% A(1,1) = x(2);
% A(1,2) = x(1);
% A(2,1) = x(3)/x(1);
% A(2,2) = x(2);
% A(2,3) = 1;
% 
% A(3,2) = A(3,2)+2*x(3);
% A(3,3) = A(3,3)+ x(2);
% A(3,4) = A(3,4)+6/x(1);

N=length(x);
A = sym(zeros(N,N));
% % For grad's system 
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
% 
% % QMBE
% 
i=N-1;
A(i,3) = (i-2)/2*(x(i-1)*~isequal((i-3)*(i-4),0))+ ...
    x(3)*x(i-3)/2*~isequal((i-5)*(i-6),0)...
    - (i+1)*(i)/2*x(i+1)/x(3);
i=N;
% A(i,2) = 0;
A(i,3) = - x(i-1)+ x(i-3)*~isequal((i-5)*(i-6),0)*x(3)/2;
A(i,4) = -3*x(i-2)*~isequal((i-4)*(i-5),0)/x(1) + 3*(i)*x(i)/x(1)/x(3);
if(N==5)
    A(i,4) = A(i,4)+ x(3);
end


A(1,1) = x(2);
A(1,2) = x(1);
A(2,1) = x(3)/x(1);
A(2,2) = x(2);
A(2,3) = 1;

A(3,2) = A(3,2)+2*x(3);
A(3,3) = A(3,3)+ x(2);
A(3,4) = A(3,4)+6/x(1);

A_l = A;
clear A_l
N = M+1;
A = sym(zeros(N,N));
% % For grad's system 
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
% 
% % QMBE
% 
i=N-1;
A(i,3) = (i-2)/2*(x(i-1)*~isequal((i-3)*(i-4),0))+x(3)*x(i-3)/2*~isequal((i-5)*(i-6),0)...
    - (i+1)*(i)/2*x(i+1)/x(3);
i=N;
% A(i,2) = 0;
A(i,3) = - x(i-1)+ x(i-3)*~isequal((i-5)*(i-6),0)*x(3)/2;
A(i,4) = -3*x(i-2)*~isequal((i-4)*(i-5),0)/x(1) + 3*(i)*x(i)/x(1)/x(3);
if(N==5)
    A(i,4) = A(i,4)+ x(3);
end


A(1,1) = x(2);
A(1,2) = x(1);
A(2,1) = x(3)/x(1);
A(2,2) = x(2);
A(2,3) = 1;

A(3,2) = A(3,2)+2*x(3);
A(3,3) = A(3,3)+ x(2);
A(3,4) = A(3,4)+6/x(1);