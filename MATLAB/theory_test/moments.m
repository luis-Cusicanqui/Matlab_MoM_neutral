function Hen = moments(mu,sigma,n)
Hen(1,1) = sym(1);
Hen(2,1) = mu;

for i=3:n+1
%    Hen(i,1) = mu*Hen(i-1) + sigma^2*(i-2)*Hen(i-2); 
   Hen(i,1) = mu*Hen(i-1) + diff(Hen(i-1),mu); 

end
simplify(Hen);
end