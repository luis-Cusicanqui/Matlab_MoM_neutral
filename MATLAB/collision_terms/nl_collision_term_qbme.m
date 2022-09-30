function rhs = nl_collision_term_qbme(u_sol,u_p,theta_p,Rcx,Ri)
%%%
%     Function used to evaluate the total collision term resulting form the
%   collision events considered the QBME model at one point in space and time. 
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
    mu = (u_p-u_sol(2))/sqrt(u_sol(3));
    sigma =sqrt(theta_p/u_sol(3));
    K  = special_integral(N,mu,sigma);
    
    rhs(1) = -Ri*u_sol(1);
    rhs(2) = u_sol(3)^(1/2)*u_sol(1) * Rcx  * K(2)/factorial(1);
    rhs(3) = u_sol(3)^(2/2)*u_sol(1) * Rcx  * K(3)/factorial(2);
    for i=4:N
        rhs(i) = u_sol(3)^((i-1)/2)*u_sol(1) * Rcx  * K(i)/factorial(i-1) - Rcx*u_sol(i) - Ri*u_sol(i);
    end
    
    rhs(N) = rhs(N) - rhs(2)*u_sol(N-1)/u_sol(1) - rhs(3)*u_sol(N-2)*~(N-2-3==0)/u_sol(1) ...
        + rhs(3)*N*u_sol(N)/(u_sol(1)*u_sol(3));
    for i=N-1:-1:5
        rhs(i) = rhs(i) - rhs(2)*u_sol(i-1)/u_sol(1) - rhs(3)*u_sol(i-2)*~(i-2-3==0)/u_sol(1);
    end
    rhs(2) = rhs(2)/u_sol(1);
    rhs(3) = rhs(3)*2/u_sol(1);
%     rhs(5) = u_sol(3)^(4/2)*u_sol(1) * Rcx  * K(5) - Rcx*u_sol(5);
end