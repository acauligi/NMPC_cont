function [MP_Prob,L_e_full,s_t] = ...
    setup_MP(n,m,...
               f,B,df,B_Jacob,state_con,u_con,...
               N,Tp,dt,...
               P,alpha,RPI_bound,...
               x_eq,u_eq,obs,J,Q,R,Name)

%% Constants

%State bounds
x_L = state_con(:,1);
x_U = state_con(:,2);

%Control bounds
u_L = u_con(:,1);
u_U = u_con(:,2);

%State cost weighting
% Q = 0*diag([0;0;1;1;1;4]);
% Q = diag(zeros(n,1));

%Control cost weighting
% R = eye(m);

%Number of collocation points
K = N;

%CGL nodes
[s_t,w] = clencurt(K); %t_t: [-1, 1] : <-> : [0, Tp]
s = fliplr(s_t); %t: [1, -1]

%% Final solution interpolation matrix

tau_full = 0:dt:Tp; 
s_e_full = (2*tau_full - Tp)/Tp; %[-1, 1]

%Lagrange polynomial evaluation at the interpolation points
% L_e = compute_Lagrange(length(s_e)-1,N,s_e,s_t);
L_e_full = compute_Lagrange(length(s_e_full)-1,N,s_e_full,s_t);

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = kron(D,eye(n));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]

% n_vars = (N+1)*(n+m);

%% Define problem

% u_eq = zeros(m,1);

x_eq_all = kron(ones(N+1,1),x_eq);
u_eq_all = kron(ones(N+1,1),u_eq);

xu_eq = [x_eq_all;u_eq_all];

Q_bar = kron(diag(w),Q); R_bar = kron(diag(w),R);
Q_tilde = Q_bar + kron(diag([zeros(N,1);(2/Tp)]),P);

F = blkdiag(Q_tilde,R_bar);
% F_pattern = sparse(F~=0);
    
%B_full = kron(eye(N+1),B);

xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L)];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U)]; 
       
MPC_cost = @(xu) (Tp/2)*(xu-xu_eq)'*F*(xu-xu_eq);% + repmat([zeros(6,1); 1/(N+1); zeros(3,1)], [N+1,1]).*xu;% + Obs_cost(xu,n,N,obs) ;
MPC_grad = @(xu) Tp*F*(xu-xu_eq);% + repmat([zeros(6,1); 1/(N+1); zeros(3,1)], [N+1 1]);% + Obs_grad(xu,n,m,N,obs);
MPC_hess = @(xu) Tp*F;

% MP_con_fnc = @(xu,Prob) MP_con(xu,Prob,n,geo_MPC.geo_Prob,N,P,D,f,B_full,Tp,obs);

c_L = [zeros(n*(N+1)+2,1); ones(obs.n_obs*(N+1),1)];
c_U = [zeros(n*(N+1),1);alpha;alpha;Inf*ones(obs.n_obs*(N+1),1)];

xu0 = zeros((n+m)*(N+1),1);

% Name = 'MP';
MP_Prob = conAssign(MPC_cost,MPC_grad,MPC_hess,[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, [],[],[],...
            'MP_con','MP_conJ',[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
        
% NMPC_Prob.SOL.optPar(10) = 1e-3;
        
MP_Prob.user.x_act = zeros(n,1);
MP_Prob.user.D = D;
MP_Prob.user.n = n;
MP_Prob.user.m = m;
MP_Prob.user.N = N;
MP_Prob.user.f = f;
MP_Prob.user.df = df;
%MP_Prob.user.B_full = B_full;
MP_Prob.user.B = B;
MP_Prob.user.B_Jacob = B_Jacob;
MP_Prob.user.Tp = Tp;
MP_Prob.user.P = P;
MP_Prob.user.x_eq = x_eq;
MP_Prob.user.obs = obs;
MP_Prob.user.J = J;

end