function A = compute_lin_A_at(x,u_p,theta_p)
%%%
% Function that evaluates the system matrix of the LIN model at a specific time on a specific cell.
%     Parameters
%     ----------
%     x: Array
%         Vector containing the moments at the specific center of the cell.
%     u_p: Double
%         Vector containing the plasma velocity at the specific center of the cell.
%     theta_p: Double
%         Vector containing the plasma temperature at the specific center of the cell.
% 
%     Returns
%     -------
%     A: numpy array
%         Matrix containing the evaluation of the system matrix of the model. Dimensions (M+1,M+1)
% 
%     written by Luis Fernando Cusicanqui Lopez
%%%
N = length(x);
A = zeros(N);
for i=2:N-1
    A(i,i) = u_p;
    A(i,i-1) = theta_p;
    A(i,i+1) = i;
end
A(1,1) = u_p;
A(1,2) = 1;
A(N,N) = u_p;
A(N,N-1) = theta_p;