function [u1,theta1] = compute_vel_and_temp(sol,u_p,theta_p)
%%%
% Computes the neutral velocity and temperature for the linear model.
%
% Parameters
% ----------
%   sol: Array
%       Matrix containing the moments of the LIN model over the whole domain.
%   u_p: Array
%       Vector containing the plasma velocity over the whole domain
%   theta_p: Array
%       Vector containing the plasma temperature over the whole domain
% 
% Returns
% -------
%   u1: Array
%       Matrix containing the neutral velocity over the whole domain
%   theta1: Array
%       Matrix containing the neutral temperature over the whole domain
% 
% Written by Luis Fernando Cusicanqui Lopez
%%%
[N,Nx,Nt] = size(sol);
    u = zeros(Nx,Nt);
    theta = zeros(Nx,Nt);

    u1 = sol(2,:,:)./sol(1,:,:) + u_p;
    theta1 = theta_p - (u1-u_p).^2 + 2.*sol(3,:,:)./sol(1,:,:);
    u(:,:) = u1(1,:,:);
    theta(:,:) = theta1(1,:,:);
end