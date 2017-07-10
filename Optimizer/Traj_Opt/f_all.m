function out = f_all(f,X,U,J,n,m,N)
    out = zeros(n*(N+1),1);
    for j = 1:N+1
        out(n*(j-1)+1:n*j) = f(X(1+n*(j-1):n*j),U(m*(j-1)+1:m*j),J);
    end
end
