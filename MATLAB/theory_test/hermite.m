function Hen = hermite(x,n)
Hen(1,1) = sym(1);
Hen(2,1) = x;

for i=3:n+1
   Hen(i,1) = x*Hen(i-1) - (i-2)*Hen(i-2); 
end
simplify(Hen);
end