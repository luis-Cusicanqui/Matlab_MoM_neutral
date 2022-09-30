function [u1,theta1] = compute_macros(rho,mom,energy)
%%%
% Computes the neutral velocity and temperature starting from the
% macroscopic quantities i.e, density momentum and energy.
%
% Parameters
% ----------
%   rho: Array
%       Vector containing the neutral density over the whole domain
%   mom: Array
%       Vector containing the neutral momentum over the whole domain
%   energy: Array
%       Vector containing the neutral energy over the whole domain
% 
% Returns
% -------
%   u1: Array
%       Vector containing the neutral velocity over the whole domain
%   theta1: Array
%       Vector containing the neutral temperature over the whole domain
% 
% Written by Luis Fernando Cusicanqui Lopez
%%%
    [Nx,Nt] = size(rho);
    u = zeros(Nx,Nt);
    theta = zeros(Nx,Nt);

    u1 = mom./rho;
    theta1 = (energy-rho.*u1.^2/2).*2./rho;
end