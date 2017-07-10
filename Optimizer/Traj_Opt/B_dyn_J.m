function f_all = B_dyn_J(df,X,U,J,m,n,N)
    %f_all = zeros(n*(N+1),m*(N+1));
    f_all = zeros(n*(N+1),(n+m)*(N+1));
    for j = 1:N+1
        %f_all(n*(j-1)+1:n*j,m*(j-1)+1:m*j) = B_Jacobian(x(n*(j-1)+1:n*j), u(m*(j-1)+1:m*j),J);
        
        B_Jacobian_local = B_Jacobian(X(n*(j-1)+1:n*j), U(m*(j-1)+1:m*j),J);
        
        f_all(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = B_Jacobian_local(:,1:n);
        f_all(n*(j-1)+1:n*j,m*(j-1)+n+1:m*j+n) = B_Jacobian_local(:,n+1:end);
    end
end