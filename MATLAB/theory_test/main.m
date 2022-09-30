clear all;
% close all;

%% PArameters
method = 'upwind';

Nt = 501;
Nx = 1001;
L = 10;
T = 4;
dx = L/(Nx-1);
dt = T/(Nt-1);
t = 0:dt:T;
x = 0:dx:L;
x = x';
N = length(x);
%%
% Initial condition
syms f(k)
f(k) = piecewise(k<1,0,k>=1 &k<=3,1,k>3,0);
initial = f(x);
U_init(1:3,:) = [initial';initial';initial'];
% X_init(:,1) = cos(x-2);
% X_init(X_init<0)= 0;
% X_init(1000:end,1) = 0;
U_sol = solve_equation(dt,dx,U_init,method,L,T,N,length(t),0,0);


%% Figures
figure(1)
hold on
plot(x,U_sol(1,:,1))
hold on
% plot(x,U_sol(1,:,200))
% plot(x,U_sol(1,:,500))
% plot(x,U_sol(1,:,1000))
% plot(x,U_sol(1,:,1500))
plot(x,U_sol(1,:,end))

title(method);

figure(2)
hold on
plot(x,U_sol(2,:,1))
hold on
% plot(x,U_sol(2,:,200))
% plot(x,U_sol(2,:,500))
% plot(x,U_sol(2,:,1000))
% plot(x,U_sol(2,:,1500))
plot(x,U_sol(2,:,end))

title(method);

figure(3)
hold on
plot(x,U_sol(3,:,1))
hold on
% plot(x,U_sol(3,:,200))
% plot(x,U_sol(3,:,500))
% plot(x,U_sol(3,:,1000))
% plot(x,U_sol(3,:,1500))
plot(x,U_sol(3,:,end))

title(method);

% 
% 
%    