function [u_aux,converged] = ...
compute_Maux(aux_prob,x_nom,u_nom,...
state_0,t_shift,Tf,Tp,N,n,m,P,f,df,dt,dt_aux)

%x_nom, u_nom: collocation point solution coefficients

%Given: collocation points for aux_problem
%Find: desired (real) time points to evaluate Nominal
tau = (1/2)*(Tf*aux_prob.s_t+Tf) + t_shift;

%Find: evaluation (scaled) time points
s_e_des = (2*tau - Tp)/Tp;

%Compute Lagrange
L_e = compute_Lagrange(length(s_e_des)-1,...
       N,s_e_des,aux_prob.s_t); %scaled nodes are the same

%Compute trajectories
x_des = zeros((N+1)*n,1);
u_des = zeros((N+1)*m,1);
for k = 1:N+1
    x_des(1+(k-1)*n:k*n) = x_nom'*L_e(:,k);
    u_des(1+(k-1)*m:k*m) = u_nom'*L_e(:,k);
end
xu_des = [x_des; u_des];

%Now we have "desired" state & control values at the nodes

%% Define Problem

Q_tilde = aux_prob.Q_bar + kron(diag([zeros(N,1);(2/Tf)]),P);
F = blkdiag(Q_tilde,aux_prob.R_bar);
aux_cost = @(xu) (Tf/2)*(xu-xu_des)'*F*(xu-xu_des);
aux_grad = @(xu) Tf*F*(xu-xu_des);

c_L = zeros(n*(N+1),1);
c_U = zeros(n*(N+1),1);

xu0 = [state_0;x_des(n+1:end);u_des];
b_L = state_0; b_U = state_0;

Name = 'Aux';
Prob = conAssign(aux_cost,aux_grad,[],[],...
        aux_prob.xu_L,aux_prob.xu_U,Name, xu0,...
            [], 0, aux_prob.Aeq,b_L,b_U,...
            'NMPC_conM','NMPC_conJM',[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
Prob.user.D = aux_prob.D;
Prob.user.n = n;
Prob.user.m = m;
Prob.user.N = N;
Prob.user.f = f;
Prob.user.df = df;
Prob.user.B_full = aux_prob.B_full;
Prob.user.Tf = Tf;

Prob = ProbCheck(Prob, 'snopt');
Result = snoptTL(Prob);

converged = Result.Inform; %GOOD: {1,2,3}

tau = 0:dt:dt_aux;
s_e = (2*tau - Tf)/Tf;
%Lagrange polynomial evaluation at the interpolation points
L_e = compute_Lagrange(length(s_e)-1,N,s_e,aux_prob.s_t);

%Compute trajectories
u_aux = zeros(size(L_e,2),m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    u_aux(:,j) = (c*L_e)';
end






end