function f_all = B_dyn_J(B_Jacobian,x,u,J,m,n,N)
    f_all = zeros(n*(N+1),m*(N+1));
    for j = 1:N+1
        f_all(n*(j-1)+1:n*j,m*(j-1)+1:m*j) = B_Jacobian(x(n*(j-1)+1:n*j), u(m*(j-1)+1:m*j),J);
    end
end