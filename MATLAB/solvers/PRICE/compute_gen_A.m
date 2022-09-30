function A_gen = compute_gen_A(u_l,u_r,u_p_left,u_p_right,theta_p_left,theta_p_right,A_comp)
%%%
% Computation of the generalized Roe matrix using three quadrature points.
%     Parameters
%     ----------
%     u_l: Array
%         Vector containing the moments at the left boundary of the current cell. Length (M+1)
%     u_r: Array
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
%     A_gen: Array
%         Resulting generalized Roe matrix. Dimensions (M+1,M+1)
% 
%     written by Luis Fernando Cusicanqui Lopez
%%%
%We assume the simple line between U_L and U_R
   s = [ 1/2-sqrt(15)/10
   1/2 
   1/2+sqrt(15)/10];
   w = [ 5/18 8/18 5/18];
   A_gen = zeros(length(u_l));
     for i=1:3
         x = u_l + s(i)*(u_r-u_l);
         u_p = u_p_left + s(i)*(u_p_right-u_p_left);
         theta_p = theta_p_left + s(i)*(theta_p_right-theta_p_left);
         A_gen =  A_gen + w(i)*A_comp(x,u_p,theta_p);
     end
end