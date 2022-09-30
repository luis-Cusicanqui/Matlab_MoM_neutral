syms v
N = size(A,1);
s = vpa(simplify(eig(A)));
p = hermite(v,N);
p = vpasolve(p(end))
s