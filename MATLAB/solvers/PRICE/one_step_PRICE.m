function x_next = one_step_PRICE(x,dx,dt,rhs,u_p,theta_p,Rcx,Ri,A_comp,lin,dirichlet)
%%%
%  Implementation of one step using PRICE method.
%     Parameters
%     ----------
%     x:  Matrix
%         Matrix containing the solution at the current timestep. Dimensions (M+1,Nx).
%     dt: Double
%         Timestep considered for simulating the model.
%     dx: Double
%         Distance between two grid points.
%     rhs: Function handle
%         Function handle that implements the collision term of the used
%         model.
%     u_p: Array
%         Vector containing the plasma velocity over the whole domain. Dimensions (Nx,1)
%     theta_p: Array
%         Vector containing the plasma temperature over the whole domain. Dimensions (Nx,1)
%     Rcx: Double
%         Double representing the charge-exchange collision rate.
%     Ri: Double
%         Double representing the ionization collision rate.
%     A_comp: Function handle
%         Function used to evaluate the system matrix of the model.
%     lin: Boolean
%         Boolean indicating whether the desired model is linear or not.
%     dirichlet: Boolean
%         Boolean indicating if boundary is dirichlet or periodic. True if dirichlet and false if periodic.
% 
%     Returns
%     -------
%     x_next: Array
%         Matrix containing the solution at the next timestep. Dimensions (M+1,Nx).
%
%     written by Luis Fernando Cusicanqui Lopez
%%%
    L = length(x);
    x_next = x;
    c = dt/dx;
%     parpool(2);
    if(dirichlet)
    if lin
      dup = zeros(L,1);
      dthetap = zeros(L,1);
      dup(1) = (u_p(2) - u_p(end))/(2*dx);
      dthetap(1) = (theta_p(2) - theta_p(end))/(2*dx);

      dup(2:L-1) = (u_p(3:end) - u_p(1:end-2))/(2*dx);
      dthetap(2:L-1) = (theta_p(3:end) - theta_p(1:end-2))/(2*dx);
      
      dup(L) = (u_p(1) - u_p(L-1))/(2*dx);
      dthetap(L) = (theta_p(1) - theta_p(L-1))/(2*dx);
      
      [A_curr_pos,A_neg_1] = compute_A_price(x(:,1),x(:,2),dx,dt,u_p(1),u_p(2),theta_p(1),theta_p(2),A_comp);
    for i=2:L-1  
        [A_next_pos,A_curr_neg] = compute_A_price(x(:,i),x(:,i+1),dx,dt,u_p(i),u_p(i+1),theta_p(i),theta_p(i+1),A_comp);
        x_next(:,i) = x(:,i) - c*(A_curr_neg*(x(:,i+1)-x(:,i)) + ...
            A_curr_pos*(x(:,i)-x(:,i-1))) ...
           +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri) ...
            +dt*lin*lin_collision(x(:,i),u_p(i),theta_p(i),dup(i),dthetap(i));%*x(:,i);
        A_curr_pos = A_next_pos;
    end  
    else
        [A_curr_pos,A_neg_1] = compute_A_price(x(:,1),x(:,2),dx,dt,u_p(1),u_p(2),theta_p(1),theta_p(2),A_comp);
        for i=2:L-1   
        [A_next_pos,A_curr_neg] = compute_A_price(x(:,i),x(:,i+1),dx,dt,u_p(i),u_p(i+1),theta_p(i),theta_p(i+1),A_comp);            
        x_next(:,i) = x(:,i) - c*(A_curr_neg*(x(:,i+1)-x(:,i)) + ...
            A_curr_pos*(x(:,i)-x(:,i-1))) ...
           +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri);
        A_curr_pos = A_next_pos;
        end
    end
    else
    if lin
      dup = zeros(L,1);
      dthetap = zeros(L,1);
      dup(1) = (u_p(2) - u_p(end))/(2*dx);
      dthetap(1) = (theta_p(2) - theta_p(end))/(2*dx);

      dup(2:L-1) = (u_p(3:end) - u_p(1:end-2))/(2*dx);
      dthetap(2:L-1) = (theta_p(3:end) - theta_p(1:end-2))/(2*dx);
      
      dup(L) = (u_p(1) - u_p(L-1))/(2*dx);
      dthetap(L) = (theta_p(1) - theta_p(L-1))/(2*dx);
      
      [A_curr_pos,A_neg_1] = compute_A_price(x(:,1),x(:,2),dx,dt,u_p(1),u_p(2),theta_p(1),theta_p(2),A_comp);
    for i=2:L-1  
        [A_next_pos,A_curr_neg] = compute_A_price(x(:,i),x(:,i+1),dx,dt,u_p(i),u_p(i+1),theta_p(i),theta_p(i+1),A_comp);
        x_next(:,i) = x(:,i) - c*(A_curr_neg*(x(:,i+1)-x(:,i)) + ...
            A_curr_pos*(x(:,i)-x(:,i-1))) ...
           +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri) ...
            +dt*lin*lin_collision(x(:,i),u_p(i),theta_p(i),dup(i),dthetap(i));%*x(:,i);
        A_curr_pos = A_next_pos;
    end  
    i=L;
    [A_next_pos,A_curr_neg] = compute_A_price(x(:,L),x(:,1),dx,dt,u_p(L),u_p(1),theta_p(L),theta_p(1),A_comp);            
    x_next(:,i) = x(:,i) - c*(A_curr_neg*(x(:,1)-x(:,i)) + ...
        A_curr_pos*(x(:,i)-x(:,i-1))) ...
       +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri) ...
       +dt*lin*lin_collision(x(:,i),u_p(i),theta_p(i),dup(i),dthetap(i));%*x(:,i);
    i=1;
    x_next(:,i) = x(:,i) - c*(A_neg_1*(x(:,i+1)-x(:,i)) + ...
        A_next_pos*(x(:,i)-x(:,end))) ...
       +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri) ...
       +dt*lin*lin_collision(x(:,i),u_p(i),theta_p(i),dup(i),dthetap(i));%*x(:,i);;
    else
        [A_curr_pos,A_neg_1] = compute_A_price(x(:,1),x(:,2),dx,dt,u_p(1),u_p(2),theta_p(1),theta_p(2),A_comp);
        for i=2:L-1   
        [A_next_pos,A_curr_neg] = compute_A_price(x(:,i),x(:,i+1),dx,dt,u_p(i),u_p(i+1),theta_p(i),theta_p(i+1),A_comp);            
        x_next(:,i) = x(:,i) - c*(A_curr_neg*(x(:,i+1)-x(:,i)) + ...
            A_curr_pos*(x(:,i)-x(:,i-1))) ...
           +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri);
        A_curr_pos = A_next_pos;
        end
        i=L;
        [A_next_pos,A_curr_neg] = compute_A_price(x(:,L),x(:,1),dx,dt,u_p(L),u_p(1),theta_p(L),theta_p(1),A_comp);            
        x_next(:,i) = x(:,i) - c*(A_curr_neg*(x(:,1)-x(:,i)) + ...
            A_curr_pos*(x(:,i)-x(:,i-1))) ...
           +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri);
        i=1;
        x_next(:,i) = x(:,i) - c*(A_neg_1*(x(:,i+1)-x(:,i)) + ...
            A_next_pos*(x(:,i)-x(:,end))) ...
           +dt*rhs(x(:,i),u_p(i),theta_p(i),Rcx,Ri);
    end
    end

end
