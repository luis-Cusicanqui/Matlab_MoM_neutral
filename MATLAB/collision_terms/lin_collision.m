function rhs = lin_collision(u_sol,u_p,theta_p,dup,dtheta)
%%%
% Function used to evaluate the extra term appearing in the collision term of the LIN model.
%     Parameters
%     ----------
%     u_sol: Array
%         Vector containing the moments at one point in space. Dimensions (M+1)
%     u_p: Double
%         Double representing the plasma velocity at the same specific position.
%     theta_p: Double
%         Double representing the plasma temperature at the same specific position.
%     dup: Double
%         Double representing the derivative of plasma velocity at the same specific position.
%     dtheta: Double
%         Double representing the derivative of plasma temperature at the same specific position.
% 
%     Returns
%     -------
%     rhs: Array
%         Vector containing the resulting extra collision term.
% 
%     written by Luis Fernando Cusicanqui Lopez
%%%
    N = length(u_sol);
    rhs = zeros(N,1);
    for i=4:N
        rhs(i) = dup*(u_p*u_sol(i-1) + theta_p*u_sol(i-2) + i*u_sol(i));
        rhs(i) = rhs(i) + dtheta*(u_p*u_sol(i-2)/2 + theta_p/2*u_sol(i-3) + i/2*u_sol(i-1));
    end
    i=1;
    rhs(1) = dup*(i*u_sol(i));
    i=2;
    rhs(2) = dup*(u_p*u_sol(i-1) + i*u_sol(i));
    rhs(2) = rhs(2) + dtheta*(i/2*u_sol(i-1));
    i=3;
    rhs(3) =dup*(u_p*u_sol(i-1) + theta_p*u_sol(i-2) + i*u_sol(i));
    rhs(3) = rhs(3) + dtheta*(u_p*u_sol(i-2)/2+ i/2*u_sol(i-1));
    rhs = -rhs;