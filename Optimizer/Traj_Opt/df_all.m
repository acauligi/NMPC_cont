function df_all = df_all(df,X,U,J,n,m,N)
    df_all = zeros(n*(N+1),(n+m)*(N+1));
    for j = 1:N+1
        % returns n*(N+1) x (n+m)*(N+1) Hessian matrix
        df_local = df(X(n*(j-1)+1:n*j),U(m*(j-1)+1:m*j),J);
        
        % first n*(N+1) columns correspond to 2nd derivatives w.r.t. X
        df_all(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = df_local(:,1:n);
        
        % last m*(N+1) columns correspond to 2nd derivatives w.r.t. U
        df_all(n*(j-1)+1:n*j,n*(N+1)+m*(j-1)+1:n*(N+1)+m*j) = df_local(:,n+1:end);
    end
end