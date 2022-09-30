clear all;
% close all;
 fprintf("Starting test: DISCONTUINITY DENSITY. \n");
  fprintf("-----------------------------------------\n\n");

% Same initial conditions but different Rcx to see the accuracy of the
% model for different.
%% PArameters
method = 'price';
place = "data/chapter5_norm/models/disc/";
test_name = place + "test_weak_disc_rho_1000_models";
type = 0;
lin = 0;
save_work = 0;


x_left = 0;
x_right = 1;
L = x_right-x_left;
T = 1e-3;

M = 10;
Ri = 0;
Rcx = 1e3;

T_ave = 1e-3;
ave=1;

%PLASMA CONSTANTS
up_amplitude = 1e0;
thetap_amplitude = 1e2;

%NEUTRAL INITIAL CONDITION CONSTANT
amplitude_rho_high =7;
amplitude_rho_low = 1;
amplitude_u_init = 1e0;
amplitude_theta_init = 1e2;

syms v
hermite_pol = hermite(v,M+1);
max_eig = up_amplitude + sqrt(thetap_amplitude)*max(vpasolve(hermite_pol(end)));

CFL = 0.7;

Nx_cells = 1000;

dx = L/(Nx_cells);

dt =1.2e-5;
t = 0:dt:T;


Nt = length(t);

Nt_ave = length(0:dt:T_ave);


x_start = x_left+dx/2;
x_end = x_right-dx/2;

Nx = floor(L./dx) + 1;

x = x_start:dx(1):x_end;

Nx = length(x);

%% Choose which simulations
simulate1 = 1;
simulate2 = 1;
simulate3 = 0;
simulate4 = 0;

sim_hme = 0;
sim_qbme = 0;
sim_lin = 1;


%% PREDEFINED Functions 
epsilon = 0.2;
syms rho(k) rho_const(k) u(k) u2(k) theta(k)
x_u_left = -0.25;
x_u_right = 0.25;
x_disc = 0.5;
rho(k) = piecewise(k<x_disc,amplitude_rho_high,k>=x_disc,amplitude_rho_low);

%% Initial conditions
U_init_hme = zeros(M+1,Nx);
U_init_hme(1,:) = double(rho(x));
U_init_hme(3,:) = amplitude_theta_init*ones(1,Nx);
U_init_hme(2,:) = amplitude_u_init*ones(1,Nx);

U_init_qbme = zeros(M+1,Nx);
U_init_qbme(1,:) = double(rho(x));
U_init_qbme(3,:) = amplitude_theta_init*ones(1,Nx);
U_init_qbme(2,:) = amplitude_u_init*ones(1,Nx);

U_init_lin = zeros(M+1,Nx);
U_init_lin(1,:) = double(rho(x));



u_p = up_amplitude*ones(1,Nx);
theta_p=thetap_amplitude*ones(1,Nx);


dirichlet = 0;
%%
if sim_hme
%     clear U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave
tic;
lin = 0;
type = 0;
[U_sol_hme,time_hme1] = solve_equation(dt,dx,CFL,U_init_hme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave,dirichlet);
[sol_hme_ave,time_hme_ave1] = solve_equation(dt,dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave,dirichlet);
[rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
toc
mom_hme_ave1 = compute_average(mom_hme,Nx,Nt_ave);
rho_hme_ave1 = compute_average(rhon_hme,Nx,Nt_ave);
energy_hme_ave1 = compute_average(energy_hme,Nx,Nt_ave);
end
if sim_qbme
tic;
lin = 0;
type = 1;
[U_sol_qbme,time_qbme1] = solve_equation(dt,dx,CFL,U_init_qbme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave,dirichlet);
[sol_qbme_ave,time_qbme_ave1] = solve_equation(dt,dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave,dirichlet);
[rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
toc
mom_qbme_ave1 = compute_average(mom_qbme,Nx,Nt_ave);
rho_qbme_ave1 = compute_average(rhon_qbme,Nx,Nt_ave);
energy_qbme_ave1 = compute_average(energy_qbme,Nx,Nt_ave);
end
if sim_lin
tic;
lin = 1;
type = 1;
[U_sol_lin,time_lin1] = solve_equation(dt,dx,CFL,U_init_lin,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave,dirichlet);
[sol_lin_ave,time_lin_ave1] = solve_equation(dt,dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave,dirichlet);
[rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
toc
mom_lin_ave1 = compute_average(mom_lin,Nx,Nt_ave);
rho_lin_ave1 = compute_average(rhon_lin,Nx,Nt_ave);
energy_lin_ave1 = compute_average(energy_lin,Nx,Nt_ave);

end

% 