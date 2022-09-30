function A = compute_nl_hme(x,u_p,theta_p)
%%%
% Function that evaluates the system matrix of the HME model at a specific time on a specific cell.
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
A = zeros(N,N);

A(1,1) = x(2);
A(1,2) = x(1);
A(2,1) = x(3)/x(1);
A(2,2) = x(2);
A(2,3) = 1;

A(3,2) = A(3,2)+2*x(3);
A(3,3) = A(3,3)+ x(2);
A(3,4) = A(3,4)+6/x(1);
for i=4:N
    A(i,i) = x(2);
    A(i,i-1) = x(3) * ~isequal((i-4),0);
    if(i<N)
        A(i,i+1) = i;
    end
    A(i,1) = -x(3)/x(1)*x(i-1)*~isequal((i-4),0);
    A(i,2) = (i)*x(i);
    A(i,3) = A(i,3)+1/2*x(3)*x(i-3)*~isequal((i-5)*(i-6),0) ...
        +(i-2)/2*x(i-1)*~isequal((i-4),0);
    A(i,4) =A(i,4) -3/x(1)*x(i-2)*~isequal((i-4)*(i-5),0);
end


A(N,2) = 0;
A(N,3) = -x(N-1) + x(N-3)*~isequal((N-5)*(N-6),0)*x(3)/2;
