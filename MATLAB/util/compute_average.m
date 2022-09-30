function u_ave = compute_average(u,Nx,Nt)
%%%
% Computes the average of a given solution over time and space.
% 
% Parameters
% ----------
% u: Array
%   Matrix containing the solution over the whole domain at every timestep
%   considered. Dimensions should be (M+1,Nx,Nt+1).
% Nx: Integer
%   Number of cells over the whole domain.
% Nt: Integer
%   Number of timesteps considered.
% Returns
% -------
% u_ave: Array
%   Vector containing the average value over time and space.
% 
% Written by Luis Fernando Cusicanqui Lopez
%%%
for i=1:Nx
    u_ave(i,:) = sum(u(i,:))/(Nt);
end
end