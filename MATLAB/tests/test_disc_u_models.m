clear all;
% close all;
 fprintf("Starting test: DISCONTUINITY DENSITY. \n");
  fprintf("-----------------------------------------\n\n");

% Experiments with same initial conditions but different Rcx to see the accuracy of the
% model for different.
%% PArameters
method = 'price';
place = "data/models/disc/";
test_name = place + "test_weak_disc_u_1000_models";

save_work = 0;

%Problem Dimensions
x_left = 0;
x_right = 1;
L = x_right-x_left;
T = 1e-2;

%Number of moments
M = 10;

%Collision terms
Ri = 0;
Rcx = [1e3 1e2 1e1 1e0];
Rcx_text = ["1e3" "1e2" "1e1" "1e0"];

T_ave = 1e-3;
ave=1;

%PLASMA CONSTANTS
up_amplitude_left = 5e0;
up_amplitude_right = 1e0; 
thetap_amplitude = 1e2;
%NEUTRAL INITIAL CONDITION CONSTANT
amplitude_rho_high =7e0;
amplitude_rho_low = 1e0;
amplitude_u_left = 5e0;
amplitude_u_right= 1e0;
amplitude_theta_init = 1e2;

syms v
hermite_pol = hermite(v,M+1);
max_eig = up_amplitude + sqrt(thetap_amplitude)*max(vpasolve(hermite_pol(end)));

%Define CFL
CFL = 0.7;

%Define number of cells
Nx_cells = 1000;

%Computation of other variables needed
dx = L/(Nx_cells);

dt = 1e-5;
t = 0:dt:T;


Nt = length(t);

Nt_ave = length(0:dt:T_ave);


x_start = x_left+dx/2;
x_end = x_right-dx/2;

Nx = floor(L./dx) + 1;

x = x_start:dx(1):x_end;

Nx = length(x);

%% Choose which simulations
simulate1 = 0;
simulate2 = 1;
simulate3 = 0;
simulate4 = 0;

sim_hme = 1;
sim_qbme = 1;
sim_lin = 1;


%% PREDEFINED Functions 
epsilon = 0.2;
syms rho(k) rho_const(k) u(k) u2(k) theta(k)
x_u_left = -0.25;
x_u_right = 0.25;
x_disc = 0.5;
rho(k) = amplitude_rho_high;
u(k) = piecewise(k<=0.5,amplitude_u_left,k>0.5,amplitude_u_right);
theta(k) = amplitude_theta_init;

syms u_equi(k) theta_equi(k)

u_equi(k) = piecewise(k<=0.5,up_amplitude_left,k>0.5,up_amplitude_right);
theta_equi(k) = thetap_amplitude;
%% Initial conditions
U_init_hme = zeros(M+1,Nx);
U_init_hme(1,:) = double(rho(x));
U_init_hme(3,:) = double(theta_equi(x));
U_init_hme(2,:) = double(u_equi(x));

U_init_qbme = zeros(M+1,Nx);
U_init_qbme(1,:) = double(rho(x));
U_init_qbme(3,:) = double(theta_equi(x));
U_init_qbme(2,:) = double(u_equi(x));

U_init_lin = zeros(M+1,Nx);
U_init_lin(1,:) = double(rho(x));


u_p = double(u_equi(x));
theta_p=double(theta_equi(x));



%%
if simulate1
    if sim_hme
    clear U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave
    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme1] = solve_equation(dt,dx,CFL,U_init_hme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type);
    [sol_hme_ave,time_hme_ave1] = solve_equation(dt,dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave1 = compute_average(mom_hme,Nx,Nt_ave);
    rho_hme_ave1 = compute_average(rhon_hme,Nx,Nt_ave);
    energy_hme_ave1 = compute_average(energy_hme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(1)+"_hme","time_hme1","time_hme_ave1",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave1','rho_hme_ave1','energy_hme_ave1');
    end
    end
    if sim_qbme
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme1] = solve_equation(dt,dx,CFL,U_init_qbme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave1] = solve_equation(dt,dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave1 = compute_average(mom_qbme,Nx,Nt_ave);
    rho_qbme_ave1 = compute_average(rhon_qbme,Nx,Nt_ave);
    energy_qbme_ave1 = compute_average(energy_qbme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(1)+"_qbme","time_qbme1","time_qbme_ave1",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave1','rho_qbme_ave1','energy_qbme_ave1');
    end
    end
    if sim_lin
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin1] = solve_equation(dt,dx,CFL,U_init_lin,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type);
    [sol_lin_ave,time_lin_ave1] = solve_equation(dt,dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(1),Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave1 = compute_average(mom_lin,Nx,Nt_ave);
    rho_lin_ave1 = compute_average(rhon_lin,Nx,Nt_ave);
    energy_lin_ave1 = compute_average(energy_lin,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(1)+"_lin",'time_lin1','time_lin_ave1','x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave1','rho_lin_ave1','energy_lin_ave1');
    end
    end
end
% 
if simulate2
    if sim_hme
    clear U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave

    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme2] = solve_equation(dt,dx,CFL,U_init_hme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(2),Ri,lin,type);
    [sol_hme_ave,time_hme_ave2] = solve_equation(dt,dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(2),Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave2 = compute_average(mom_hme,Nx,Nt_ave);
    rho_hme_ave2 = compute_average(rhon_hme,Nx,Nt_ave);
    energy_hme_ave2 = compute_average(energy_hme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(2)+"_hme","time_hme2","time_hme_ave2",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave2','rho_hme_ave2','energy_hme_ave2');
    end
    end
    if sim_qbme
    
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme2] = solve_equation(dt,dx,CFL,U_init_qbme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(2),Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave2] = solve_equation(dt,dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(2),Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave2 = compute_average(mom_qbme,Nx,Nt_ave);
    rho_qbme_ave2 = compute_average(rhon_qbme,Nx,Nt_ave);
    energy_qbme_ave2 = compute_average(energy_qbme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(2)+"_qbme","time_qbme2","time_qbme_ave2",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave2','rho_qbme_ave2','energy_qbme_ave2');
    end
    end
    if sim_lin
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin2] = solve_equation(dt,dx,CFL,U_init_lin,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(2),Ri,lin,type);
    [sol_lin_ave,time_lin_ave2] = solve_equation(dt,dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(2),Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave2 = compute_average(mom_lin,Nx,Nt_ave);
    rho_lin_ave2 = compute_average(rhon_lin,Nx,Nt_ave);
    energy_lin_ave2 = compute_average(energy_lin,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(2)+"_lin",'time_lin2','time_lin_ave2','x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave2','rho_lin_ave2','energy_lin_ave2');
    end
    end
end
if simulate3
    if sim_hme
    clear U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave

    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme3] = solve_equation(dt,dx,CFL,U_init_hme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(3),Ri,lin,type);
    [sol_hme_ave,time_hme_ave3] = solve_equation(dt,dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(3),Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave3 = compute_average(mom_hme,Nx,Nt_ave);
    rho_hme_ave3 = compute_average(rhon_hme,Nx,Nt_ave);
    energy_hme_ave3 = compute_average(energy_hme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(3)+"_hme","time_hme3","time_hme_ave3",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave3','rho_hme_ave3','energy_hme_ave3');
    end
    end
    
    if sim_qbme
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme3] = solve_equation(dt,dx,CFL,U_init_qbme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(3),Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave3] = solve_equation(dt,dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(3),Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave3 = compute_average(mom_qbme,Nx,Nt_ave);
    rho_qbme_ave3 = compute_average(rhon_qbme,Nx,Nt_ave);
    energy_qbme_ave3 = compute_average(energy_qbme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(3)+"_qbme","time_qbme3","time_qbme_ave3",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave3','rho_qbme_ave3','energy_qbme_ave3');
    end
    end
    if sim_lin
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin3] = solve_equation(dt,dx,CFL,U_init_lin,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(3),Ri,lin,type);
    [sol_lin_ave,time_lin_ave3] = solve_equation(dt,dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(3),Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave3 = compute_average(mom_lin,Nx,Nt_ave);
    rho_lin_ave3 = compute_average(rhon_lin,Nx,Nt_ave);
    energy_lin_ave3 = compute_average(energy_lin,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(3)+"_lin",'time_lin3','time_lin_ave3','x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave3','rho_lin_ave3','energy_lin_ave3');
    end
    end
end

if simulate4
    if sim_hme
    clear U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave    
    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme4] = solve_equation(dt,dx,CFL,U_init_hme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(4),Ri,lin,type);
    [sol_hme_ave,time_hme_ave4] = solve_equation(dt,dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(4),Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave4 = compute_average(mom_hme,Nx,Nt_ave);
    rho_hme_ave4 = compute_average(rhon_hme,Nx,Nt_ave);
    energy_hme_ave4 = compute_average(energy_hme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(4)+"_hme","time_hme4","time_hme_ave4",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave4','rho_hme_ave4','energy_hme_ave4');
    end
    end
    if sim_qbme
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme4] = solve_equation(dt,dx,CFL,U_init_qbme,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(4),Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave4] = solve_equation(dt,dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(4),Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave4 = compute_average(mom_qbme,Nx,Nt_ave);
    rho_qbme_ave4 = compute_average(rhon_qbme,Nx,Nt_ave);
    energy_qbme_ave4 = compute_average(energy_qbme,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(4)+"_qbme","time_qbme4","time_qbme_ave4",'x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave4','rho_qbme_ave4','energy_qbme_ave4');
    end
    end
    if sim_lin
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin4] = solve_equation(dt,dx,CFL,U_init_lin,method,4,T,Nx,Nt,-2,0,u_p,theta_p,Rcx(4),Ri,lin,type);
    [sol_lin_ave,time_lin_ave4] = solve_equation(dt,dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave,-2,0,u_p,theta_p,Rcx(4),Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave4 = compute_average(mom_lin,Nx,Nt_ave);
    rho_lin_ave4 = compute_average(rhon_lin,Nx,Nt_ave);
    energy_lin_ave4 = compute_average(energy_lin,Nx,Nt_ave);
    if save_work
    save(test_name+"_rcx_eq_"+Rcx_text(4)+"_lin",'time_lin4','time_lin_ave4','x','Rcx','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave4','rho_lin_ave4','energy_lin_ave4');
    end
    end
end

  fprintf("FINISHED TEST: SMOOTH RHO MODELS\n\n");
fprintf("-----------------------------------------\n\n");
