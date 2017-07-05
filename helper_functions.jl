using JuMP

#   ChebyshevMatrices.m: form cosine-spaced vector of independent variable, 
#                           Chebyshev differentiation and integration matrices
#
#   Ross Allen, ASL, Stanford University
#   Michael Colonno, ADL, Stanford University 
#
#   Started:        6/27/11
#   Last Updated:   Feb 13, 2014
#
#   Inputs:         N               order of polynomial approximation 
#                   t               reverse(original) time points 1 ---> -1 
#
#   Outputs:        D               differentiation matrix (non-dim)
#                   err             nonzero if an error occurs
#
#
#   References: 
#              Q. Gong, I. M. Ross, F. Fahroo, "A Chebyshev Pseudospectral
#               Method for Nonlinear Constrained Optimal Control Problems",
#               Joint 48th IEEE Conference on Decision and Control and 28th
#               Chinese Control Conference, Shanghai, P.R. China, December
#               16-18, 2009
#
#              L. Trfethen, "Spectral Methods in MATLAB", SIAM 2000,
#               Page 128
#
#   Note:
#               Francisco Capistran found problems/inaccuracies in using
#               the integration matrix. Use with caution.
#################################################################

function ChebyshevDiffMatrix(N,t)
  # differentiation matrix
  D = zeros(Float64, N+1,N+1);
  
  #compute using t (reverse ordered)
  for k in 0:N
      i_k = k+1;
      c_k = k == 0 || k==N ? 2 : 1;
      
      for j in 0:N
          i_j = j+1;
          c_j = j == 0 || j==N ? 2 : 1;
      
          if (j == k)     
              D[i_k,i_j] = -0.5*t[i_k]/(1 - t[i_k]^2);       
          else
              D[i_k,i_j] = c_k/c_j*((-1)^(j+k))/(t[i_k] - t[i_j]);
          end
      end
  end
  
  # fix corners
  D[1,1] = (2*(N^2) +1)/6; 
  D[N+1,N+1] = -D[1,1];
  
  # adjust for forward ordered time
  D = -D; 
  
  return D
end
#################################################################





# CLENCURT   nodes x (Chebyshev points) and weights w 
#            for Clenshaw-Curtis quadrature
#
#   Reference: L. Trfethen, "Spectral Methods in MATLAB", SIAM 2000,
#              Page 128

function clencurt(K)
  theta = pi*(0:K)/K;
  x = cos(theta)[end:-1:1]; #-1 ---> 1 (t)
  
  if K%2 == 0
    w = zeros(div(K,2)+1,1); #s:0 ---> N/2
    w[1] = 1/(K^2-1); #w_0
    
    w_pre = [0.5+0.5/(1-K^2)*cos(pi*n) for n in 1:K/2]; # s:1--->N/2

    [w[s+1] = 4/K*w_pre[s] + 4/K*sum([1/(1-4*j^2)*cos(2*pi*j*s/K) for j in 1:div(K,2)-1]) for s in 1:div(K,2)];

    w_net = [w; w[end-1:-1:1]]
  end 

  return x,w_net
end




function compute_Lagrange(K,N,t,t_nodes)
	L_e = ones(N+1,K+1); #up to order N

	#[L_e[k+1,:] = L_e[k+1,:].*((t-t_nodes[m+1])/(t_nodes[k+1]-t_nodes[m+1])) for m in 0:N for k in 0:N if m!=k];

	for m in 0:N
		for k in 0:N
			if m!=k
				#L_e[k+1,:] = L_e[k+1,:].*((t-t_nodes[m+1])/(t_nodes[k+1]-t_nodes[m+1]));
				for j in 1:length(t)
					L_e[k+1,j] = L_e[k+1,j]*((t[j]-t_nodes[m+1])/(t_nodes[k+1]-t_nodes[m+1]));
				end
			end
		end
	end

	return L_e
end




function df_all(state, N_collocation, A, B, n, m)
  out = [];
  for i in 1:N_collocation+1
    x = state[n*(i-1)+1:n*i];
    u = state[n*(N+1)+m*(i-1)+1:n*(N+1)+m*i];

    out = [out; A*x + B*u]
  end
  return out
end
