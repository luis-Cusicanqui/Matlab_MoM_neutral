function x_next = one_step_FV(x,dx,dt,rhs)

L = length(x);

x_next = x;
A = 1;
C = A*dt/dx;
for i=2:L-1
    x_next(:,i) = x(:,i) - dt/dx*(compute_flux_at(x(:,i))- compute_flux_at(x(:,i-1)));  
end
x_next(:,1) = x(:,1) - dt/dx*(compute_flux_at(x(:,1))- compute_flux_at(x(:,end)));  
end
