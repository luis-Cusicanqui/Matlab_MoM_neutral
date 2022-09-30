function rhs = neutral_collision_term(u_sol,u_p,theta_p,Kn,Ri)
%%%
%     Function used to evaluate the collision term in mono-atomic gases. 
%     
%     Parameters
%     ----------
%     u_sol: Array
%         Vector containing the moments at one point in space. Dimensions (M+1)
%     u_p: Double
%         Double representing the plasma velocity at the same specific position.
%     theta_p: Double
%         Double representing the plasma temperature at the same specific position.
%     Kn: Double
%         Knudsen number considerd
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
    rhs(1) = 0;
    Rcx = u_sol(1)/Kn;
    for i=4:N
        rhs(i) = -Rcx*u_sol(i);
    end
end