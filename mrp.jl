using JuMP, Ipopt, PyPlot

include("helper_functions.jl");

mod = Model(solver=IpoptSolver(print_level=0));

pdot(p,w) = [0.25*(p[1]^2-p[2]^2-p[3]^2+1)*w[1] +  0.5*(p[1]*p[2]-p[3])*w[2] + 0.5*(p[2]-p[1]*p[3])*w[3]; 0.5*(p[1]*p[2]+p[3])*w[1] 0.25*(p[2]^2 - p[1]^2-p[3]^2)*w[2] + 0.5*(p[2]*p[3]-p[1])*w[3]; 0.5*(p[1]*p[3]-p[2])*w[1] + 0.5*(p[1]+p[2]*p[3])*w[2] + 0.25*(p[3]^2-p[2]^2-p[1]^2+1)*w[3]];
pdot(p,w) = [p[1]*w[1]; p[2]*w[2]; p[3]*w[3]];

# Sim parameters
mass = 6.0;
w_max = 1.2;

N = 120;
n = 3;
m = 3;

dt = 0.01;
Tp = 5;

# BCs
x0 = [-1.0; 1.0; 1.0];
xf = [0;0;0];
u_eq = zeros(m,1);

# Constraints
state_constr_lo = -Inf*ones(n,1);
state_constr_hi = -state_constr_lo;

ctrl_constr_lo = -w_max*ones(m,1);
ctrl_constr_hi = -state_constr_lo;

Xmin = [kron(ones(N+1), state_constr_lo); kron(ones(N+1), ctrl_constr_lo)];
Xmax = [kron(ones(N+1), state_constr_hi); kron(ones(N+1), ctrl_constr_hi)];

Q = zeros(n,n);     R = eye(m);

# CGL nodes
s_t, w = clencurt(N);
s = s_t[end:-1:1];

# Final solution interpolation matrix
tau_full = 0:dt:Tp;
s_e_full = [(2*tau - Tp)/Tp for tau in tau_full];

# Lagrange polynomical evaluation at the interpolation pts
Le = compute_Lagrange(length(s_e_full)-1, N, s_e_full, s_t);

# Get differentiation matrix
D = ChebyshevDiffMatrix(N,s);
D = kron(ChebyshevDiffMatrix(N,s), eye(n));

@variable(mod, Xmin[i] <= state[i=1:(n+m)*(N+1)] <= Xmax[i]);

# LQR cost
#Q_bar = sparse(kron(diagm([w[i] for i in 1:length(w)]), Q));
#R_bar = sparse(kron(diagm([w[i] for i in 1:length(w)]), R));
#Q_tilde = blkdiag(Q_bar, R_bar);
#lqr_cost = state'*Q_tilde*state;

lqr_cost = 0;
for i = 1:N+1
  x = state[n*(i-1)+1:n*i];
  u = state[n*(N+1)+m*(i-1)+1:n*(N+1)+m*i];
  lqr_cost += ((x-xf)'*diagm(w[i]*ones(n))*Q*(x-xf) + (u-u_eq)'*diagm(w[i]*ones(m))*R*(u-u_eq));
end

# Objective: minimize quadratic LQR cost of trajectory
@objective(mod, Min, lqr_cost[1]);

# BC constraints 
@constraints(mod, begin
    state[1:n] .== x0;
    state[n*N+1:n*(N+1)] .== xf;
end)

# Dynamics constraint
JuMP.register(mod, :pdot, 2, pdot, autodiff=true);
JuMP.register(mod, :df_all, 5, df_all, autodiff=true);
#@NLconstraint(mod, df_all(pdot, state, N, n, m) - 2/Tp*D*state[1:n*(N+1)] == zeros(n*(N+1)))

## Solve problem
#status = solve(mod);
#println("Solver status: ", status)
#
## Unpack solution
#soln = getvalue(state);
#state_out = zeros(size(Le,2),n);
#ctrl_out = zeros(size(Le,2),m);
#
#if status == :Optimal
#  for i = 1:n
#    c = soln[i:n:n*(N+1) - (n-i)];
#  	state_out[:,i] = [(c'*Le[:,i])[1] for i in 1:size(Le,2)];
#  end
#  
#  for i = 1:m
#    c = soln[n*(N+1)+i:m:end-(m-i)];
#  	ctrl_out[:,i] = [(c'*Le[:,i])[1] for i in 1:size(Le,2)];
#  end
#end
#
#state_out = state_out';
#ctrl_out = ctrl_out';
#
#plot(state_out[1,:], state_out[2,:]);