clear all
close all

syms x
N=4;
hermitefunf = hermite(x,N);
hermite1 = hermite(x,10);

weight(x) = 1/(sqrt(2*pi))*exp(-x.^2/2);
f(x) = hermite1(9);
point = vpasolve(hermitefunf(end));
w = weight(point);
w2 = compute_weights(N,hermitefunf(end-1),point);

sum = 0;
for i=1:N
   sum = sum + w2(i)*f(point(i)); 
end


function w = compute_weights(N,hermiteN_1,points)
    syms x
    h2 = hermiteN_1*hermiteN_1;
    w = factorial(N-1)/N./subs(h2,x,points);
end