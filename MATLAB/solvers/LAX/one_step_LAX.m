function x_next = one_step_LAX(x,dx,dt,rhs)

    L = length(x);

    x_next = x;
    C = dt/dx;

    for i=2:L-1
        x_next(i) = x(i) - C/2*(compute_flux_at(x(i+1)) - compute_flux_at(x(i-1))) ...
            + C^2/2*(compute_A_at(0.5*(x(i) + x(i+1)))*(compute_flux_at(x(i+1)) -compute_flux_at(x(i))) ...
            - compute_A_at(0.5*(x(i-1) + x(i)))*(compute_flux_at(x(i)) -compute_flux_at(x(i-1))));
    end
    i=1;
    x_next(1) =x(i) - C/2*(compute_flux_at(x(i+1)) - compute_flux_at(x(end))) ...
            + C^2/2*(compute_A_at(0.5*(x(i) + x(i+1)))*(compute_flux_at(x(i+1)) -compute_flux_at(x(i))) ...
            - compute_A_at(0.5*(x(end) + x(i)))*(compute_flux_at(x(i)) -compute_flux_at(x(end)))); 
    i=L;
    x_next(L) =   x(i) - C/2*(compute_flux_at(x(1)) - compute_flux_at(x(i-1))) ...
            + C^2/2*(compute_A_at(0.5*(x(i) + x(1)))*(compute_flux_at(x(1)) -compute_flux_at(x(i))) ...
            - compute_A_at(0.5*(x(i-1) + x(i)))*(compute_flux_at(x(i)) -compute_flux_at(x(i-1))));
end