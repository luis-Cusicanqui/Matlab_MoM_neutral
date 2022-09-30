function rhs = lin_collision_term(u_sol,u_p,theta_p,Rcx,Ri)
%%%
%     Function used to evaluate the total collision term resulting form the
%   collision events considered the LIN model at one point in space and time. 
%     
%     Parameters
%     ----------
%     u_sol: Array
%         Vector containing the moments at one point in space. Dimensions (M+1)
%     u_p: Double
%         Double representing the plasma velocity at the same specific position.
%     theta_p: Double
%         Double representing the plasma temperature at the same specific position.
%     Rcx: Double
%         Collision rate considered.
%     Ri: Double
%         Collision rate considered.
% 
%     Returns
%     -------
%     rhs: Array
%         Vector containing the resulting ionization term.
% 
%     written by Luis Fernando Cusicanqui Lopez
%%%
    N = length(u_sol);
    rhs = zeros(N,1);
    rhs(1) = -Ri*u_sol(1);
    rhs(2) = - Rcx*u_sol(2) - Ri*u_sol(2);
    rhs(3) = - Rcx*u_sol(3) - Ri*u_sol(3);
    for i=4:N
        rhs(i) = - Rcx*u_sol(i) - Ri*u_sol(i);
    end
end