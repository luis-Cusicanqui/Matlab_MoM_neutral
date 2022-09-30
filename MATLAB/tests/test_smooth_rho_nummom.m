clear all;
close all;
 fprintf("Starting test: SMOOTH DIFFERENT MOMENTS. \n");
fprintf("-----------------------------------------\n\n");

%% PArameters
method = 'price';
savefig = 0;
plot_fig = 0;
save_work = 1;
place = "data/moments/smooth/";
test_name = place + "test_smooth_rho_1000_rcx_1e0_mom";

x_left = 0;
x_right = 1;
L = x_right-x_left;
T = 1e-2;

M = [4 5 6 7 8];
Ri = 0;
Rcx = 1e0;

T_ave = 1e-3;
ave=1;

%PLASMA CONSTANTS
up_amplitude = 1e0;
thetap_amplitude = 1e2;

%NEUTRAL INITIAL CONDITION CONSTANT
amplitude_rho_ave =7e0;
amplitude_rho_amp = 3e0;
amplitude_u_init = 1e0;
amplitude_theta_init = 1e2;

syms v
hermite_pol = hermite(v,M(5)+1);
max_eig(1) = double(up_amplitude + sqrt(thetap_amplitude)*max(vpasolve(hermite_pol(M(1)+2))));
max_eig(2) = double(up_amplitude + sqrt(thetap_amplitude)*max(vpasolve(hermite_pol(M(2)+2))));
max_eig(3) = double(up_amplitude + sqrt(thetap_amplitude)*max(vpasolve(hermite_pol(M(3)+2))));
max_eig(4) = double(up_amplitude + sqrt(thetap_amplitude)*max(vpasolve(hermite_pol(M(4)+2))));
max_eig(5) = double(up_amplitude + sqrt(thetap_amplitude)*max(vpasolve(hermite_pol(M(5)+2))));

CFL = 0.7;


Nx_cells = 1000;

dx = L/(Nx_cells);

%dt = dx.*CFL./double(max_eig);
dt = 1.25e-5* ones(1,5);
%dt = 5e-6* ones(1,5);
CFL = dt./dx.*double(max_eig);
t1 = 0:dt(1):T;
t2 = 0:dt(2):T;
t3 = 0:dt(3):T;
t4 = 0:dt(4):T;
t5 = 0:dt(5):T;

Nt = [length(t1);length(t2);length(t3);length(t4);length(t5)];

Nt_ave = [length(0:dt(1):T_ave);length(0:dt(2):T_ave); ...
    length(0:dt(3):T_ave);length(0:dt(4):T_ave);length(0:dt(5):T_ave)];


x_start = x_left+dx/2;
x_end = x_right-dx/2;

Nx = floor(L./dx) + 1;

x = x_start:dx(1):x_end;

Nx = length(x);

clear t1 t2 t3 t4 t5

%% Choose which simulations
simulate1 = 0;
simulate2 = 0;
simulate3 = 0;
simulate4 = 0;
simulate5 = 1;

sim_hme = 1;
sim_qbme = 1;
sim_lin = 1;




%% PREDEFINED Functions 
syms rho(k) rho_const(k) u(k) u2(k) theta(k)
x_u_left = -0.25;
x_u_right = 0.25;
x_disc = 0.5;
rho(k) = amplitude_rho_ave + amplitude_rho_amp*sin(2*pi*k);
%% Initial condition
u_p = up_amplitude*ones(1,Nx);
theta_p=thetap_amplitude*ones(1,Nx);


%%
if simulate1
    if sim_hme
    U_init_hme = zeros(M(1)+1,Nx);
    U_init_hme(1,:) = double(rho(x));
    U_init_hme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_hme(2,:) = amplitude_u_init*ones(1,Nx);


    
    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme1] = solve_equation(dt(1),dx,CFL,U_init_hme,method,4,T,Nx,Nt(1),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_hme_ave,time_hme_ave1] = solve_equation(dt(1),dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave(1),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave1 = compute_average(mom_hme,Nx,Nt_ave(1));
    rho_hme_ave1 = compute_average(rhon_hme,Nx,Nt_ave(1));
    energy_hme_ave1 = compute_average(energy_hme,Nx,Nt_ave(1));
    if save_work
    save(test_name+"_eq_"+M(1)+"_hme",'time_hme1','time_hme_ave1','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave1','rho_hme_ave1','energy_hme_ave1');
    end
    end
    if sim_qbme
    U_init_qbme = zeros(M(1)+1,Nx);
    U_init_qbme(1,:) = double(rho(x));
    U_init_qbme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_qbme(2,:) = amplitude_u_init*ones(1,Nx);
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme1] = solve_equation(dt(1),dx,CFL,U_init_qbme,method,4,T,Nx,Nt(1),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave1] = solve_equation(dt(1),dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave(1),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave1 = compute_average(mom_qbme,Nx,Nt_ave(1));
    rho_qbme_ave1 = compute_average(rhon_qbme,Nx,Nt_ave(1));
    energy_qbme_ave1 = compute_average(energy_qbme,Nx,Nt_ave(1));
    if save_work
    save(test_name+"_eq_"+M(1)+"_qbme",'time_qbme1','time_qbme_ave1','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave1','rho_qbme_ave1','energy_qbme_ave1');
    end

    end

    if sim_lin
    U_init_lin = zeros(M(1)+1,Nx);
    U_init_lin(1,:) = double(rho(x));
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin1] = solve_equation(dt(1),dx,CFL,U_init_lin,method,4,T,Nx,Nt(1),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_lin_ave,time_lin_ave1] = solve_equation(dt(1),dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave(1),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave1 = compute_average(mom_lin,Nx,Nt_ave(1));
    rho_lin_ave1 = compute_average(rhon_lin,Nx,Nt_ave(1));
    energy_lin_ave1 = compute_average(energy_lin,Nx,Nt_ave(1));
    if save_work
    save(test_name+"_eq_"+M(1)+"_lin",'time_lin1','time_lin_ave1','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave1','rho_lin_ave1','energy_lin_ave1');
    end
    end
    clear U_init_hme U_init_lin U_init_qbme U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave
end
% 
if simulate2
    if sim_hme
    U_init_hme = zeros(M(2)+1,Nx);
    U_init_hme(1,:) = double(rho(x));
    U_init_hme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_hme(2,:) = amplitude_u_init*ones(1,Nx);


    
    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme2] = solve_equation(dt(2),dx,CFL,U_init_hme,method,4,T,Nx,Nt(2),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_hme_ave,time_hme_ave2] = solve_equation(dt(2),dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave(2),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave2 = compute_average(mom_hme,Nx,Nt_ave(2));
    rho_hme_ave2 = compute_average(rhon_hme,Nx,Nt_ave(2));
    energy_hme_ave2 = compute_average(energy_hme,Nx,Nt_ave(2));
    if save_work
    save(test_name+"_eq_"+M(2)+"_hme",'time_hme2','time_hme_ave2','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave2','rho_hme_ave2','energy_hme_ave2');
    end
    end
    if sim_qbme
    U_init_qbme = zeros(M(2)+1,Nx);
    U_init_qbme(1,:) = double(rho(x));
    U_init_qbme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_qbme(2,:) = amplitude_u_init*ones(1,Nx);

    
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme2] = solve_equation(dt(2),dx,CFL,U_init_qbme,method,4,T,Nx,Nt(2),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave2] = solve_equation(dt(2),dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave(2),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave2 = compute_average(mom_qbme,Nx,Nt_ave(2));
    rho_qbme_ave2 = compute_average(rhon_qbme,Nx,Nt_ave(2));
    energy_qbme_ave2 = compute_average(energy_qbme,Nx,Nt_ave(2));
    if save_work
    save(test_name+"_eq_"+M(2)+"_qbme",'time_qbme2','time_qbme_ave2','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave2','rho_qbme_ave2','energy_qbme_ave2');
    end
    end
    if sim_lin
    U_init_lin = zeros(M(2)+1,Nx);
    U_init_lin(1,:) = double(rho(x));
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin2] = solve_equation(dt(2),dx,CFL,U_init_lin,method,4,T,Nx,Nt(2),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_lin_ave,time_lin_ave2] = solve_equation(dt(2),dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave(2),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave2 = compute_average(mom_lin,Nx,Nt_ave(2));
    rho_lin_ave2 = compute_average(rhon_lin,Nx,Nt_ave(2));
    energy_lin_ave2 = compute_average(energy_lin,Nx,Nt_ave(2));
    if save_work
    save(test_name+"_eq_"+M(2)+"_lin",'time_lin2','time_lin_ave2','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave2','rho_lin_ave2','energy_lin_ave2');
    end
    end
    end
    clear U_init_hme U_init_lin U_init_qbme U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave

end
if simulate3
    if sim_hme
    U_init_hme = zeros(M(3)+1,Nx);
    U_init_hme(1,:) = double(rho(x));
    U_init_hme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_hme(2,:) = amplitude_u_init*ones(1,Nx);

    

    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme3] = solve_equation(dt(3),dx,CFL,U_init_hme,method,4,T,Nx,Nt(3),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_hme_ave,time_hme_ave3] = solve_equation(dt(3),dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave(3),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave3 = compute_average(mom_hme,Nx,Nt_ave(3));
    rho_hme_ave3 = compute_average(rhon_hme,Nx,Nt_ave(3));
    energy_hme_ave3 = compute_average(energy_hme,Nx,Nt_ave(3));
    if save_work
    save(test_name+"_eq_"+M(3)+"_hme",'time_hme3','time_hme_ave3','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave3','rho_hme_ave3','energy_hme_ave3');
    end
    end
    if sim_qbme
    U_init_qbme = zeros(M(3)+1,Nx);
    U_init_qbme(1,:) = double(rho(x));
    U_init_qbme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_qbme(2,:) = amplitude_u_init*ones(1,Nx);

    
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme3] = solve_equation(dt(3),dx,CFL,U_init_qbme,method,4,T,Nx,Nt(3),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave3] = solve_equation(dt(3),dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave(3),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave3 = compute_average(mom_qbme,Nx,Nt_ave(3));
    rho_qbme_ave3 = compute_average(rhon_qbme,Nx,Nt_ave(3));
    energy_qbme_ave3 = compute_average(energy_qbme,Nx,Nt_ave(3));
    if save_work
    save(test_name+"_eq_"+M(3)+"_qbme",'time_qbme3','time_qbme_ave3','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave3','rho_qbme_ave3','energy_qbme_ave3');
    end
    end
    if sim_lin
    U_init_lin = zeros(M(3)+1,Nx);
    U_init_lin(1,:) = double(rho(x));
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin3] = solve_equation(dt(3),dx,CFL,U_init_lin,method,4,T,Nx,Nt(3),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_lin_ave,time_lin_ave3] = solve_equation(dt(3),dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave(3),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave3 = compute_average(mom_lin,Nx,Nt_ave(3));
    rho_lin_ave3 = compute_average(rhon_lin,Nx,Nt_ave(3));
    energy_lin_ave3 = compute_average(energy_lin,Nx,Nt_ave(3));
    if save_work
    save(test_name+"_eq_"+M(3)+"_lin",'time_lin3','time_lin_ave3','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave3','rho_lin_ave3','energy_lin_ave3');
    end
    end
    clear U_init_hme U_init_lin U_init_qbme U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave

end

if simulate4
    if sim_hme
    U_init_hme = zeros(M(4)+1,Nx);
    U_init_hme(1,:) = double(rho(x));
    U_init_hme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_hme(2,:) = amplitude_u_init*ones(1,Nx);

    

    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme4] = solve_equation(dt(4),dx,CFL,U_init_hme,method,4,T,Nx,Nt(4),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_hme_ave,time_hme_ave4] = solve_equation(dt(4),dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave(4),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave4 = compute_average(mom_hme,Nx,Nt_ave(4));
    rho_hme_ave4 = compute_average(rhon_hme,Nx,Nt_ave(4));
    energy_hme_ave4 = compute_average(energy_hme,Nx,Nt_ave(4));
    if save_work
    save(test_name+"_eq_"+M(4)+"_hme",'time_hme4','time_hme_ave4','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave4','rho_hme_ave4','energy_hme_ave4');
    end
    end
    if sim_qbme
    U_init_qbme = zeros(M(4)+1,Nx);
    U_init_qbme(1,:) = double(rho(x));
    U_init_qbme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_qbme(2,:) = amplitude_u_init*ones(1,Nx);

    
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme4] = solve_equation(dt(4),dx,CFL,U_init_qbme,method,4,T,Nx,Nt(4),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave4] = solve_equation(dt(4),dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave(4),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave4 = compute_average(mom_qbme,Nx,Nt_ave(4));
    rho_qbme_ave4 = compute_average(rhon_qbme,Nx,Nt_ave(4));
    energy_qbme_ave4 = compute_average(energy_qbme,Nx,Nt_ave(4));
    if save_work
    save(test_name+"_eq_"+M(4)+"_qbme",'time_qbme4','time_qbme_ave4','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave4','rho_qbme_ave4','energy_qbme_ave4');
    end
    end
    if sim_lin
    U_init_lin = zeros(M(4)+1,Nx);
    U_init_lin(1,:) = double(rho(x));
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin4] = solve_equation(dt(4),dx,CFL,U_init_lin,method,4,T,Nx,Nt(4),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_lin_ave,time_lin_ave4] = solve_equation(dt(4),dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave(4),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave4 = compute_average(mom_lin,Nx,Nt_ave(4));
    rho_lin_ave4 = compute_average(rhon_lin,Nx,Nt_ave(4));
    energy_lin_ave4 = compute_average(energy_lin,Nx,Nt_ave(4));
    if save_work
    save(test_name+"_eq_"+M(4)+"_lin",'time_lin4','time_lin_ave4','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave4','rho_lin_ave4','energy_lin_ave4');
    end
    end
    end
    clear U_init_hme U_init_lin U_init_qbme U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave
    
end

if simulate5
    if sim_hme
    U_init_hme = zeros(M(5)+1,Nx);
    U_init_hme(1,:) = double(rho(x));
    U_init_hme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_hme(2,:) = amplitude_u_init*ones(1,Nx);

    

    tic;
    lin = 0;
    type = 0;
    [U_sol_hme,time_hme5] = solve_equation(dt(5),dx,CFL,U_init_hme,method,4,T,Nx,Nt(5),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_hme_ave,time_hme_ave5] = solve_equation(dt(5),dx,CFL,U_sol_hme(:,:,end),method,4,T_ave,Nx,Nt_ave(5),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_hme,mom_hme,energy_hme] = compute_quantities(sol_hme_ave,lin,u_p,theta_p);
    toc
    mom_hme_ave5 = compute_average(mom_hme,Nx,Nt_ave(5));
    rho_hme_ave5 = compute_average(rhon_hme,Nx,Nt_ave(5));
    energy_hme_ave5 = compute_average(energy_hme,Nx,Nt_ave(5));
    if save_work
    save(test_name+"_eq_"+M(5)+"_hme",'time_hme5','time_hme_ave5','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_hme_ave5','rho_hme_ave5','energy_hme_ave5');
    end
    end
    if sim_qbme
    U_init_qbme = zeros(M(5)+1,Nx);
    U_init_qbme(1,:) = double(rho(x));
    U_init_qbme(3,:) = amplitude_theta_init*ones(1,Nx);
    U_init_qbme(2,:) = amplitude_u_init*ones(1,Nx);

    
    tic;
    lin = 0;
    type = 1;
    [U_sol_qbme,time_qbme5] = solve_equation(dt(5),dx,CFL,U_init_qbme,method,4,T,Nx,Nt(5),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_qbme_ave,time_qbme_ave5] = solve_equation(dt(5),dx,CFL,U_sol_qbme(:,:,end),method,4,T_ave,Nx,Nt_ave(5),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_qbme,mom_qbme,energy_qbme] = compute_quantities(sol_qbme_ave,lin,u_p,theta_p);
    toc
    mom_qbme_ave5 = compute_average(mom_qbme,Nx,Nt_ave(5));
    rho_qbme_ave5 = compute_average(rhon_qbme,Nx,Nt_ave(5));
    energy_qbme_ave5 = compute_average(energy_qbme,Nx,Nt_ave(5));
    if save_work
    save(test_name+"_eq_"+M(5)+"_qbme",'time_qbme5','time_qbme_ave5','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_qbme_ave5','rho_qbme_ave5','energy_qbme_ave5');
    end
    end
    if sim_lin
    U_init_lin = zeros(M(5)+1,Nx);
    U_init_lin(1,:) = double(rho(x));
    tic;
    lin = 1;
    type = 1;
    [U_sol_lin,time_lin5] = solve_equation(dt(5),dx,CFL,U_init_lin,method,4,T,Nx,Nt(5),-2,0,u_p,theta_p,Rcx,Ri,lin,type);
    [sol_lin_ave,time_lin_ave5] = solve_equation(dt(5),dx,CFL,U_sol_lin(:,:,end),method,4,T_ave,Nx,Nt_ave(5),-2,0,u_p,theta_p,Rcx,Ri,lin,type,ave);
    [rhon_lin,mom_lin,energy_lin] = compute_quantities(sol_lin_ave,lin,u_p,theta_p);
    toc
    mom_lin_ave5 = compute_average(mom_lin,Nx,Nt_ave(5));
    rho_lin_ave5 = compute_average(rhon_lin,Nx,Nt_ave(5));
    energy_lin_ave5 = compute_average(energy_lin,Nx,Nt_ave(5));
    if save_work
    save(test_name+"_eq_"+M(5)+"_lin",'time_lin5','time_lin_ave5','x','Rcx','T_ave','M','CFL','T','theta_p','u_p','dt','dx','mom_lin_ave5','rho_lin_ave5','energy_lin_ave5');
    end
    end
    clear U_init_hme U_init_lin U_init_qbme U_sol_hme U_sol_qbme U_sol_lin sol_hme_ave sol_qbme_ave sol_lin_ave
    
end

fprintf("ENDING test: SMOOTH DIFFERENT MOMENTS. \n");
fprintf("-----------------------------------------\n\n");



