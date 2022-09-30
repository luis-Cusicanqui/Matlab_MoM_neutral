function Gn = special_integral(N,mu,sigma)
%%%
% Function that implements the recursion used for the computation of the integrals of the following form:
%         int(He_i(x)NormalDist(mu,sigma^2))dx
%     for i=0..N.
%     Parameters
%     ---------
%     N: Int
%         Highest polynomial degree considered
%     mu: Double
%         Mean value of the normal distribution
%     sigma: Double
%         Standard variation of the normal distribution
% 
%     Returns
%     -------
%     Gn: Array
%         Vector containing on the entry the result of the integral.
% 
%     written by Luis Fernando Cusicanqui Lopez
%%%
    Gn(1) = 1;
    Gn(2) = mu;
    for i=3:N
        Gn(i) = mu*Gn(i-1) + (i-2)*Gn(i-2)*(sigma^2-1);
    end
end