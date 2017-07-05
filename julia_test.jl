using JuMP, Ipopt, PyPlot

include("helper_functions.jl");

mod = Model(solver=IpoptSolver(print_level=0));

# Sim parameters
mass = 6.0;
F_max = 2.5;

N = 120;
n = 6;
m = 3;

dt = 0.01;
Tp = 28;

# BCs
x0 = [-4.4;-5;0;1.3;0;0];

x_eq = [4.5;4.5;0;0;0;0];
u_eq = zeros(m,1);

# Dynamics and constraints
A = [zeros(3,3) eye(3); zeros(3,6)];
B = [zeros(3,3); eye(3)];

state_constr_lo = [-10*ones(3,1); -5*ones(3,1)];
state_constr_hi = -state_constr_lo;

ctrl_constr_lo = -F_max*ones(m,1);
ctrl_constr_hi = -state_constr_lo;

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

@variables(mod, begin
    state[1:(n+m)*(N+1)]
end)

# LQR cost
lqr_cost = 0;
for i = 1:N+1
  x = state[n*(i-1)+1:n*i];
  u = state[n*(N+1)+m*(i-1)+1:n*(N+1)+m*i];
  lqr_cost += ((x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq));
end

# Objective: minimize quadratic LQR cost of trajectory
@objective(mod, Min, lqr_cost[1]);

# BC constraints 
@constraints(mod, begin
    state[1:n] .== x0;
    state[n*N+1:n*(N+1)] .== x_eq;
end)

# Dynamics constraint
@constraint(mod, df_all(state, N, A, B, n, m) - 2/Tp*D*state[1:n*(N+1)] .== zeros(n*(N+1)))

# Solve problem
status = solve(mod);
println("Solver status: ", status)

# Unpack solution
soln = getvalue(state);
state_out = zeros(size(Le,2),n);
ctrl_out = zeros(size(Le,2),m);

if status == :Optimal
  for i = 1:n
    c = soln[i:n:n*(N+1) - (n-i)];
  	state_out[:,i] = [(c'*Le[:,i])[1] for i in 1:size(Le,2)];
  end
  
  for i = 1:m
    c = soln[n*(N+1)+i:m:end-(m-i)];
  	ctrl_out[:,i] = [(c'*Le[:,i])[1] for i in 1:size(Le,2)];
  end
end

state_out = state_out';
ctrl_out = ctrl_out';
