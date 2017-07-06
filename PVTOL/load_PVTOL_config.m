%% Constants mass = 0.486;
mass = 6.0;
w_max = 0.5;
F_max = 2.5;
T_max = 2.5;

%% Obstacle info

obs_loc = [[3;-4],...
           [0.7;-3],...
           [-1;-0.5],...
           [2.5;-0.5],...
           [-4;1],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2.5;-4.5],...
           [-2;4]];
obs_loc_mpc = [[0.7;-3],...
           [-1;-0.5],...
           [2.5;-0.5],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2.5;-4.5]];
obs_rad = [1,0.9,0.8,1.2,1,0.9,0.5,1,0.6];
obs_rad_mpc = [0.9,0.8,1.2,0.9,0.5,1];
obs_rad = [];
obs = struct('n_obs',length(obs_rad),'pos',obs_loc,'r',obs_rad);
obs_mpc = struct('n_obs',length(obs_rad_mpc),'pos',obs_loc_mpc,'r',obs_rad_mpc);

%% Setup Metric
% load 'metric_PVTOL_vectorized.mat';
% 
% % W_fnc = @(x) W_mat(x);
% % dW_fnc = @(x) {dW_p_mat(x),dW_vy_mat(x)};
% 
% W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly_fnc);
% dW_fnc = @(x) {dw_poly_p_fnc(x), dw_poly_vy_fnc(x)};
% 
% % sigma_ThBw = 0.3296;
% sigma_ThBw = 0.3185;
% lambda =  0.8283;
% ctrl_bound = 6.00;
% n_W = [3,4];

%% x = [x y th xd yd thd]'
% n = 6; m = 3;
% J = 0.00383;    g = 9.81;   
len = 0.25;
% f  = @(x) [x(4)*cos(x(3)) - x(5)*sin(x(3));
%            x(4)*sin(x(3)) + x(5)*cos(x(3));
%            x(6);
%            x(6)*x(5)-g*sin(x(3));
%            -x(6)*x(4)-g*cos(x(3));
%            0];
%        
% B = [zeros(1,4),1/mass, len/J;
%      zeros(1,4),1/mass,-len/J]';
% 
% df = @(x) [0,0,-x(4)*sin(x(3))-x(5)*cos(x(3)),cos(x(3)),-sin(x(3)),0;
%            0,0, x(4)*cos(x(3))-x(5)*sin(x(3)),sin(x(3)), cos(x(3)),0;
%            zeros(1,5),1;
%            0,0,-g*cos(x(3)), 0, x(6), x(5);
%            0,0, g*sin(x(3)), -x(6), 0, -x(4);
%            zeros(1,6)];
% 
% B_w = [zeros(1,3),1,0,0;
%        zeros(1,3),0,1,0]';
%    
% state_constr_low = -[5.5;5.5;pi/4;2;1;pi/3];%+euc_bound;
% state_constr = [state_constr_low, -state_constr_low];
% ctrl_constr = 10*[-ones(m,1) ones(m,1)];
%            
% x_eq = [4.5;4.5;0;0;0;0];
% u_eq = [0.5*mass*g; 0.5*mass*g; 0.5*mass*g]; 
% 
% test_state = [-4.4;-5;0;1.3;0;0];

%% x = [rx ry rz vx vy vz]
n = 6;  m = 3;
J = diag([5,5,10]);

f  = @(x,J) [x(4:6);
          zeros(3,1)];
       
df  = @(x,J) [zeros(3) eye(3);
                zeros(3,6)];
   
%B = @(u) 1/mass*[zeros(3,3); eye(3)]*u;
B = @(x,u,J) [zeros(3,1); u(1);u(2);u(3)];
B_Jacobian = @(x,u,J) [zeros(3,3); eye(3)];
%B_Jacobian = @(x,u,J) [zeros(3,m+n); zeros(3,n) eye(m)];

state_constr_low = -[5.5;5.5;pi/4;2;1;pi/3];
state_constr = [state_constr_low, -state_constr_low];
ctrl_constr = F_max*[-ones(m,1) ones(m,1)];
           
x_eq = [4.5;4.5;0;0;0;0];
u_eq = zeros(m,1);

test_state = [-4.4;-5;0;1.3;0;0];
filename = 'double_integrator_mat.mat';

%% x = [rx ry rz vx vy vz r]
% n = 7;  m = 3;
% J = diag([5,5,10]);
% 
% f  = @(x,J) [x(7)*x(4); x(7)*x(5); x(7)*x(6);...
%           zeros(4,1)];
%        
% df  = @(x,J) [zeros(3) x(7)*eye(3) x(4:6);
%                 zeros(3,7)];
%    
% %B = @(u) 1/mass*[zeros(3,3); eye(3)]*u;
% B = @(x,u,J) [zeros(3,1); u(1);u(2);u(3); 0];
% B_Jacobian = @(x,u,J) [zeros(3,3); eye(3); zeros(1,3)];
% B_Jacobian = @(x,u,J) [zeros(3,n+m); zeros(3,n) eye(m) zeros(3,1); zeros(1,n+m)];
% 
% state_constr_low = -[5.5;5.5;pi/4;2;1;pi/3;0];
% state_constr_hi = [5.5;5.5;pi/4;2;1;pi/3;30];
% state_constr = [state_constr_low, state_constr_hi];
% ctrl_constr = F_max*[-ones(m,1) ones(m,1)];
%            
% x_eq = [4.5;4.5;0;0;0;0;0];
% u_eq = zeros(m,1);
% 
% test_state = [-4.4;-5;0;1.3;0;0;4];
% filename = 'double_integrator_fft_mat.mat';

%% x = [p1 p2 p3 wx wy wz]'     u = [Tx Ty Tz]'
% J = diag([5,5,10]);
% n = 6;  m = 3;
% 
% f  = @(x,J) [0.25*x(4)*(x(1)^2-x(2)^2-x(3)^2+1) + 0.5*x(5)*(x(1)*x(2)-x(3)) + 0.5*x(6)*(x(2)+x(1)*x(3));
%             0.5*x(4)*(x(3)+x(1)*x(2)) + 0.25*x(5)*(x(2)^2-x(1)^2-x(3)^2+1) + 0.5*x(6)*(x(2)*x(3)-x(1));
%             0.5*x(4)*(x(1)*x(3)-x(2)) + 0.5*x(5)*(x(1)+x(2)*x(3)) + 0.25*x(6)*(x(3)^2-x(2)^2-x(1)^2+1);
%            (J(2,2) - J(3,3))*x(5)*x(6)/J(1,1);
%            (J(3,3) - J(1,1))*x(4)*x(6)/J(2,2);
%            (J(1,1) - J(2,2))*x(4)*x(5)/J(3,3)];
%        
% 
% df  = @(x,J) [0.5*(x(1)*x(4)+x(2)*x(5)+x(3)*x(6)),0.5*(x(6)-x(2)*x(4)+x(1)*x(5)),0.5*(x(1)*x(6)-x(3)*x(4)-x(5)),0.25*(x(1)^2-x(2)^2-x(3)^2+1),0.5*(x(1)*x(2)-x(3)),0.5*(x(2)+x(1)*x(3));
%             0.5*(x(2)*x(4)-x(1)*x(5)-x(5)),0.5*(x(2)*x(5)+x(1)*x(4)+x(3)*x(5)),0.5*(x(4)-x(3)*x(5)+x(2)*x(6)),0.5*(x(3)+x(1)*x(2)),0.25*(x(2)^2-x(1)^2-x(3)^2+1),0.5*(x(2)*x(3)-x(1));
%             0.5*(x(5)-x(1)*x(6)+x(3)*x(4)),0.5*(x(3)*x(5)-x(2)*x(6)-x(4)),0.5*(x(3)*x(6)+x(1)*x(4)+x(2)*x(5)),0.5*(x(1)*x(3)-x(2)),0.5*(x(1)+x(2)*x(3)),0.25*(x(3)^2-x(1)^2-x(2)^2+1);
%            0,0,0,0, (J(2,2) - J(3,3))*x(6)/J(1,1), (J(2,2) - J(3,3))*x(5)/J(1,1);
%            0,0,0, (J(3,3) - J(1,1))*x(6)/J(2,2), 0, (J(3,3) - J(1,1))*x(4)/J(2,2);
%            0,0,0, (J(1,1) - J(2,2))*x(5)/J(3,3), (J(1,1) - J(2,2))*x(4)/J(3,3), 0];
% state_constr_low = [-1*ones(3,1); -w_max*ones(3,1)];
% state_constr = [state_constr_low, -state_constr_low];
% ctrl_constr = 10*[-ones(m,1) ones(m,1)];
%            
% x_eq = zeros(n,1);
% u_eq = zeros(m,1);
% 
% test_state = [quat2mrp([-0.5  0.5 0.5 -0.5])'; 0.05*ones(3,1)];
% filename = 'mrp_T_input_mat.mat';

%% x = [p1 p2 p3]'    u = [wx wy wz]'
% n = 3;  m = 3;
% J = diag([5,5,10]);
% 
% f  = @(x) zeros(3,1);
% df  = @(x) zeros(3,3);
% 
% B = @(x,u,J) [0.25*(x(1)^2-x(2)^2-x(3)^2+1), 0.5*(x(1)*x(2)-x(3)), 0.5*(x(2)-x(1)*x(3));...
%             0.5*(x(1)*x(2)+x(3)), 0.25*(x(2)^2 - x(1)^2-x(3)^2), 0.5*(x(2)*x(3)-x(1));...
%             0.5*(x(1)*x(3)-x(2)), 0.5*(x(1)+x(2)*x(3)), 0.25*(x(3)^2-x(2)^2-x(1)^2+1)]*u;
% 
% B_Jacobian = @(x,u,J) [0.25*(x(1)^2-x(2)^2-x(3)^2+1), 0.5*(x(1)*x(2)-x(3)), 0.5*(x(2)-x(1)*x(3));...
%             0.5*(x(1)*x(2)+x(3)), 0.25*(x(2)^2 - x(1)^2-x(3)^2), 0.5*(x(2)*x(3)-x(1));...
%             0.5*(x(1)*x(3)-x(2)), 0.5*(x(1)+x(2)*x(3)), 0.25*(x(3)^2-x(2)^2-x(1)^2+1)];
% 
% th_max = 358*pi/180;    th_tange = tan(th_max/4);
% state_constr_low = -th_tange*ones(n,1);
% state_constr = [state_constr_low, -state_constr_low];
% ctrl_constr = w_max*[-ones(m,1) ones(m,1)];
%            
% x_eq = zeros(n,1);
% u_eq = zeros(m,1);
% 
% test_state = quat2mrp([-0.5  0.5 0.5 -0.5])';
% filename = 'mrp_w_input_mat.mat';

%% Dynamics and cost
f_true = f;
B_true = B;
B_w_true = B;

Q = zeros(n); R = eye(m);
   
%% Bounds

w_max = 0.1;

% M_ccm = W_upper\eye(n);
% d_bar = (w_max*sigma_ThBw/lambda);
d_bar = 0.0385;
% ctrl_bound = ctrl_bound*w_max;
% euc_bound = d_bar*sqrt(diag(W_upper));
% 
% In = eye(n);
% M_ccm_pos_unscaled = ((In(1:2,:)*W_upper*In(1:2,:)')\eye(2));
% M_ccm_pos = (1/d_bar^2)*((In(1:2,:)*W_upper*In(1:2,:)')\eye(2));
% [U_pos,S_pos,V_pos] = svd(M_ccm_pos);
%     
% %Rescale ellipsoids by obstacle + robot radius
% M_obs = zeros(2,2,obs.n_obs);
% for i = 1:obs.n_obs
%     S_new = (sqrt(S_pos\eye(2)) + (obs_rad(i)+len)*eye(2))^2\eye(2);
%     M_obs(:,:,i) = U_pos*S_new*V_pos';
% end
% obs.M_obs = M_obs;
% 
% M_obs_mpc = zeros(2,2,obs_mpc.n_obs);
% for i = 1:obs_mpc.n_obs
%     S_new = (sqrt(S_pos\eye(2)) + (obs_mpc.r(i)+len)*eye(2))^2\eye(2);
%     M_obs_mpc(:,:,i) = U_pos*S_new*V_pos';
% end
% obs_mpc.M_obs = M_obs_mpc;

P = 2.5*eye(n);
alpha = 1e-3;