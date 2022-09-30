function A = hermite_gen(n)
    A = zeros(n+1);
    A(1,1) = 1;
    A(2,2) = 1;
    for i=3:n+1
       A(i,1) = -(i-2)*A(i-2,1);
       for j=2:i
        A(i,j) = A(i-1,j-1) - (i-2)*A(i-2,j);
       end
    end
end