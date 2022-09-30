function [A_pos,A_neg] = compute_A_price(x_left,x_right,dx,dt,u_p_left,u_p_right,theta_p_left,theta_p_right,A_comp)
%%%
% Computation of the generalized Roe matrix using three quadrature points.
%     Parameters
%     ----------
%     x_left: Array
%         Vector containing the moments at the left boundary of the current cell. Length (M+1)
%     x_right: Array
%         Vector containing the moments at the right boundary of the current cell. Length (M+1)
%     u_p_left: Double
%         Variable containing the plasma velocity at center of the left cell of the current cell.
%     u_p_right: Double
%         Variable containing the plasma velocity at center of the right cell of the current cell.
%     theta_p_left: Double
%         Variable containing the plasma temperature at center of the left cell of the current cell.
%     theta_p_right: Double
%         Variable containing the plasma temperature at center of the right cell of the current cell.
%     A_comp: Function handle
%         Function used to evaluate the system matrix of the model.
%        
%     Returns
%     -------
%     A_pos: Array
%         Matrix containing the A- in the PRICE scheme
%     A_neg: Array
%          Matrix containing the A- in the PRICE scheme
%
%     written by Luis Fernando Cusicanqui Lopez
%%%    
c1 = dx/dt;
    c2 = dt/dx;
    A_temp = compute_gen_A(x_left,x_right,u_p_left,u_p_right,theta_p_left,theta_p_right,A_comp);
    N = length(x_left);
    Im = eye(N);
       A_pos = 2*A_temp + Im*c1 + c2*(A_temp)^2; 
       A_neg = 2*A_temp - Im*c1 - c2*(A_temp)^2; 
    A_pos = 1/4*A_pos;
    A_neg = 1/4*A_neg;
end