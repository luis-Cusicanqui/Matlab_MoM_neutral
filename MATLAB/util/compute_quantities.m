function [rho,momentum,energy] = compute_quantities(sol,linear,u_p,theta_p)
%%%
% Computes the macroscopic quantities, i.e, density, momentum and energy; 
% starting from the moments and the
% plasma quantities.
% 
% Parameters
% ----------
%   sol: Array
%       Matrix containing the moments of the LIN model over the whole domain.
%   u_p: Array
%       Vector containing the plasma velocity over the whole domain
%   theta_p: Array
%       Vector containing the plasma temperature over the whole domain
%   linear: boolean
%       Variable indicating whether the model used is linear or non-linear.
% Returns
% -------
%   rho: Array
%       Matrix containing the neutral density over the whole domain
%   momentum: Array
%       Matrix containing the neutral momentum over the whole domain
%   energy: Array
%       Matrix containing the neutral energy over the whole domain
% 
% Written by Luis Fernando Cusicanqui Lopez
%%%
   [N,Nx,Nt] = size(sol);
   rho = zeros(Nx,Nt);
   momentum = zeros(Nx,Nt);
   energy = zeros(Nx,Nt);
   u_n = zeros(Nx,Nt);
   theta_n = zeros(Nx,Nt);
   rho(:,:) = sol(1,:,:);
   if linear
        [u_n(:,:),theta_n(:,:)] = compute_vel_and_temp(sol,u_p,theta_p);
   else
        u_n(:,:) = sol(2,:,:);
        theta_n(:,:) = sol(3,:,:);
   end
   momentum(:,:) = rho.*u_n;
   energy(:,:) = rho.*u_n.^2/2 +rho.*theta_n/2; 
end