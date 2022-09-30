function [u_sol_fin,time_needed] = solve_equation(dt,dx,CFL,u_init,method,L,T,Nx,Nt,x_left,t_init,u_p,theta_p,Rcx,Ri,linear,type,ave,boundary)
%%%
% Driver function that sets the chooses parameters depending on the inputs. It calls the time-stepping method to solve
%     the moment model.
%     PARAMETERS
%     ----------
%     dt: Double
%         Timestep considered for simulating the model.
%     dx: Double
%         Distance between two grid points.
%     CFL: Double
%         CFL number used.
%     u_init: Array
%         Array containing the initial conditions of the problem. Dimensions should be (Nx,M+1)
%     method: Str
%         String indicating the method.      
%     L: Double
%         Length of the considered domain 
%     T: Double
%         Total simulation time.
%     Nx: Int
%         Number of grid points considered.
%     Nt: Int
%         Number of time instances considered.
%     x_left: Double
%         x-coordinate of the left boundary of domain.
%     t_init: Double
%         Initial time.
%     u_p: Array
%         Vector containing the plasma velocity over the whole domain. Dimensions (Nx,1)
%     theta_p: Array
%         Vector containing the plasma temperature over the whole domain. Dimensions (Nx,1)
%     Rcx: Double
%         Double represeting the charge-exchange collision rate.
%     Ri: Double
%         Double representing the ionization collision rate.
%     linear: Boolean
%         Boolean indicating whether the desired model is linear or not.
%     type: Int
%         If linear=False, then type is a integer that chooses between HME and QBME model.
%         0 if HME and 1 if QBME both for plasma edge simulations
%         2 if HME and 3 if QBME both for rarefied gases simulations
%     ave: Boolean
%         Boolean indicating whether all the time instances should be stored or
%         not.
%     boundary: Boolean
%         Boolean indicating if boundary is dirichlet or periodic. True if dirichlet and false if periodic.
%       
%     Returns
%     -------
%     u_sol_fin: Array
%         Array designated to store the solution at different time instances. IF all time instances are stored
%       then the dimensions should be (M+1,Nx,Nt+1)
%     time_needed: Double
%         Variable containing the total computation time of the simulation.
%     written by Luis Fernando Cusicanqui Lopez
%%%
N = size(u_init,1);
if(nargin<17)
    error('Too few arguments given')
end
if(nargin==17)
    ave=0;
end
if(nargin<=18)
    boundary=0;
end
%% Parameters
if linear
   A_comp =  @compute_lin_A_at;
    rhs = @lin_collision_term;
    lin = 1;   
    dirichlet = boundary;
    info = "linear";
else
switch (type)
    case 0
    A_comp = @compute_nl_hme;
    rhs = @nl_collision_term_hme;
    lin = 0;
    dirichlet = boundary;
    info = "hme";  
    case 1
    A_comp = @compute_nl_qbme;
    rhs = @nl_collision_term_qbme;
    lin = 0;  
    dirichlet = boundary;
    info = "qbme";
    case 2
    A_comp = @compute_nl_hme;
    rhs = @neutral_collision_term;
    lin = 0;
    dirichlet = boundary;
    info = "hme";
    case 3
    A_comp = @compute_nl_qbme;
    rhs = @neutral_collision_term;
    lin = 0;
    dirichlet = boundary;
    info = "qbme";
end
end
%% Method
switch(method)
    case 'lax'
        method = @one_step_LAX;
    case 'upwind'
        method = @one_step_FV;
    case 'price'
        method = @one_step_PRICE;
end

%% Main For-loop
fprintf("SIMULATION QUANTITIES\n");
fprintf("----------------------\n\n");
fprintf("Model used:"+info+"\n")
fprintf("CFL used: %d\n",CFL)
fprintf("Number of cells:%d\n",Nx);
fprintf("Number of instances: %d\n",Nt);
fprintf("Number of moments: %d \n",N-1);
fprintf("Collision rate Rcx: %d \n",Rcx)
fprintf("Collision rate Ri: %d \n",Ri)
fprintf("----------------------\n");

  [u_sol_fin,time_needed] = foward_euler(method,u_init,Nx,Nt,dx,dt,rhs,u_p,theta_p,Rcx,Ri,A_comp,lin,ave,dirichlet);

    fprintf("Time required for simulation: %d\n",time_needed);
    fprintf("SIMULATION FINISHED \n\n");
    fprintf("----------------------\n\n");

