function conJ = MP_conJ(xu,Prob)%,D,n,m,N,df,B_full,Tp,P)
%Dynamics, and terminal
n = Prob.user.n;
N = Prob.user.N;
m = Prob.user.m;
J = Prob.user.J;

no = Prob.user.obs.n_obs;

global US_A;
global geodesic_MPC;

conJ = zeros(n*(N+1)+2+no*(N+1),(n+m)*(N+1));

%% Dynamics constraints

% conJ(1:n*(N+1),1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D - ...
%             df_all(Prob.user.df,xu(1:n*(N+1)),J,n,N);
% 
% conJ(1:n*(N+1), :) = conJ(1:n*(N+1), :)-B_dyn_J(Prob.user.B_Jacob,xu(1:n*(N+1)),xu(n*(N+1)+1:end),J,m,n,N);
X = xu(1:n*(N+1));
U = xu(n*(N+1)+1:end);
conJ(1:n*(N+1), :) = [(2/Prob.user.Tp)*Prob.user.D zeros(n*(N+1),m*(N+1))]-...
                        df_all(Prob.user.df,X,U,J,n,m,N);


%conJ(1:n*(N+1), n*(N+1)+1:end) = -B_dyn_J(Prob.user.B_Jacob,xu(1:n*(N+1)),xu(n*(N+1)+1:end),J,m,n,N);
%conJ(1:n*(N+1),n*(N+1)+1:end) = -Prob.user.B_full;

%% Initial RPI constraint

% conJ(end-1,1:n) = NaN;

% w_poly = geodesic_MPC.W.w_poly_fnc(xu(1:n));
% M = (geodesic_MPC.W.W_eval(w_poly))\eye(n);
% conJ(n*(N+1)+1,1:n) = -2*US_A'*M;
conJ(n*(N+1)+1,1:n) = 2*(Prob.user.P*(xu(1:n)-Prob.user.x_act))';


%% Terminal constraint
conJ(n*(N+1)+2,n*N+1:n*(N+1)) = 2*(Prob.user.P*(xu(n*N+1:n*(N+1))-Prob.user.x_eq))';

%% Obstacle constraints

% cJ_obs = zeros(no*(N+1),(n+m)*(N+1));

for i = 1:no
    o_pos = Prob.user.obs.pos(:,i);
    Mo = Prob.user.obs.M_obs(:,:,i);
    for k = 1:N+1
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
%         cJ_obs((i-1)*(N+1)+k,1+(k-1)*n:2+(k-1)*n) = -2*(o_pos-x_k)'*Mo;
        conJ(n*(N+1)+2+(i-1)*(N+1)+k,1+(k-1)*n:2+(k-1)*n) = -2*(o_pos-x_k)'*Mo;
    end
end

% conJ(n*(N+1)+3:end,:) = cJ_obs;

% conJ = conJ;


end

