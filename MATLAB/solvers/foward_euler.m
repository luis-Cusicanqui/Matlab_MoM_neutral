function [u_sol_fin,time_needed] = foward_euler(method,u_init,Nx,Nt,dx,dt,rhs,u_p,theta_p,Rcx,Ri,A_comp,lin,ave,dirichlet)
%%%
%   Implementation of a foward-euler time-stepping method used to simulate the moment model.
%     Parameters
%     ----------
%     method: Function handle
%         Function handle that computes the simulation scheme. For example
%         the PRICE scheme can be used.
%     u_init: Array
%         Array containing the initial conditions of the problem. Dimensions should be (Nx,M+1)
%     Nx: Int
%         Number of grid points considered.
%     Nt: Int
%         Number of time instances considered. 
%     dt: Double
%         Timestep considered for simulating the model.
%     dx: Double
%         Distance between two grid points.
%     rhs: Function handle
%         Function handle that implements the collision term of the used
%         model.
%     u_p: Array
%         Vector containing the plasma velocity over the whole domain. Dimensions (Nx,1)
%     theta_p: Array
%         Vector containing the plasma temperature over the whole domain. Dimensions (Nx,1)
%     Rcx: Double
%         Double representing the charge-exchange collision rate.
%     Ri: Double
%         Double representing the ionization collision rate.
%     A_comp: Function handle
%         Function used to evaluate the system matrix of the model.
%     lin: Boolean
%         Boolean indicating whether the desired model is linear or not.
%     ave: Boolean
%         Boolean indicating if boundary is dirichlet or periodic. True if dirichlet and false if periodic.
%     dirichlet: Boolean
%         Boolean indicating if boundary is dirichlet or periodic. True if dirichlet and false if periodic.
%       
%     Returns
%     -------
%     u_sol_fin: Array
%         Array designated to store the solution at different time instances. IF all timesteps are stored
%       then the dimensions should be (M+1,Nx,Nt+1)
%     time_needed: Double
%         Variable containing the total computation time of the simulation.
% 
%     written by Luis Fernando Cusicanqui Lopez
%%%
    N = size(u_init,1);
    u_sol_fin(:,:,1) = u_init; 
    u_sol = zeros(N,Nx,2); %columns is x at one specific time t
    u_sol(:,:,1) = u_init;
    tic
    for iter=1:Nt-1
       if(mod(iter,2)==0)
            index_next = 1;
            index_curr = 2;
       else
            index_next = 2;
            index_curr = 1;
       end
    u_sol(:,:,index_next) = method(u_sol(:,:,index_curr),dx,dt,rhs,u_p,theta_p,Rcx,Ri,A_comp,lin,dirichlet); 
    if ave
        u_sol_fin(:,:,end+1) = u_sol(:,:,index_next);
    elseif mod(iter,100)==0 || iter==Nt-1
        u_sol_fin(:,:,end+1) = u_sol(:,:,index_next);
    end
    
    end
    time_needed = toc;
end