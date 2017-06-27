function f_all = B_dyn(B,x,u,J,m,n,N)
    f_all = zeros(n*(N+1),1);
    for j = 1:N+1
        f_all(n*(j-1)+1:n*j) = B(x(n*(j-1)+1:n*j),u(m*(j-1)+1:m*j),J);
    end
end