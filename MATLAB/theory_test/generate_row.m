function row_nmin1 = generate_row(n,coeff_n)
    row_nmin1 = zeros(1,length(coeff_n));
    for i=1:n
        row_nmin1(i) = coeff_n(i+1) *(i)/(n);
    end