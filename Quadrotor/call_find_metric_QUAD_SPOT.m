clear all; close all; clc;
yalmip('clear');

warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Constants

n = 9;
% lambda = 1;
g = 9.81;
ccm_eps = 1e-1;

r_lim = pi/6;
p_lim = pi/6;
vx_lim = 1;
vy_lim = 1;
vz_lim = 1;
th_lim_low = 0.5*g;
th_lim_high = 2*g;

%% Define problem

% lambda_range = linspace(0.1,4.0,13);
% lambda_range = (1/100)*round(lambda_range*100);
% lambda_range = 0.75;
% euc_bounds = NaN(length(lambda_range),1);
% d_bars = NaN(length(lambda_range),1);
% cond_bound = NaN(length(lambda_range),1);
% 
% eps = 0.1;
% condn_prev = 17;
% return_metric = 0;
% 
% for ll = 1:length(lambda_range)
%     lambda = lambda_range(ll);
%     
%     fprintf('**********\n');
%     fprintf('lambda: %f\n', lambda);
%     solved = 0;
%     
%     %Determine upper bound
%     cond_l = condn_prev;
%     cond_u = 1.2*condn_prev;
%     while (~solved)
%         fprintf(' cond_u: %.2f: ', cond_u);
%         [sos_prob,~,~] = find_metric_QUAD_SPOT_89(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,...
%             cond_u,lambda,ccm_eps,return_metric);
%         if (sos_prob == 0)
%             solved = 1;
%             fprintf('feasible \n');
%         else
%             %shift up condition number range
%             fprintf('\n');
%             cond_l = cond_u;
%             cond_u = 1.2*cond_u;
%         end
%     end
%     if (solved)
%         euc_bounds(ll) = sqrt(cond_u)/lambda;
%         fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
%     else
%         continue;
%     end
%     
%     %Now do bisection search
%     while(cond_u - cond_l >= eps)
%         condn = (cond_l+cond_u)/2;
%         fprintf(' cond: %.2f', condn);
%         try
%             [sos_prob,w_lower,w_upper] = find_metric_QUAD_SPOT_89(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,...
%             condn,lambda,ccm_eps,return_metric);
%         catch
%             sos_prob = 1; 
%         end
%         if (sos_prob == 0)
%             fprintf(' feasible\n');
%             
%             euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;
%             d_bars(ll) = sqrt(double(1/w_lower))/lambda;
%             
%             cond_u = condn;
%         else
%             fprintf(' infeasible\n');
%             cond_l = condn;
%         end
%     end
%     condn_prev = cond_u;
%     cond_bound(ll) = cond_u;
%     disp('Euc_bound:'); disp(euc_bounds(ll));
%     disp('d_bar:'); disp(d_bars(ll));
%     fprintf('**********\n');
%     
% end
% 
% pause;

%% Pick a solution

% load metric_FLR_comp_final.mat;
lambda = 0.75;
condn = 21;
return_metric = 1;

[sos_prob,w_lower,w_upper] = find_metric_QUAD_SPOT_89(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,...
            condn,lambda,ccm_eps,return_metric);

%% Compute aux control bound
load 'metric_QUAD_vectorized.mat';
disp('Checking CCM conditions and Computing control bound...');
% lambda = 0.98*lambda;

Bw = [zeros(3);
    eye(3);
    zeros(3)];
B_perp = [eye(6);
    zeros(3,6)];

ctrl_N = 15;
r_range = linspace(-r_lim, r_lim, ctrl_N);
p_range = linspace(-p_lim, p_lim, ctrl_N);
th_range = linspace(th_lim_low, th_lim_high, 2*ctrl_N);

b_T = @(x) [sin(x(8)); -cos(x(8))*sin(x(7)); cos(x(8))*cos(x(7))];

f_mat = @(x) [[x(4);x(5);x(6)];
    [0;0;g] - b_T(x)*x(9);
    zeros(3,1)];

%gradients
db_T_q =@(x) [0, cos(x(8));
    -cos(x(7))*cos(x(8)), sin(x(7))*sin(x(8));
    -sin(x(7))*cos(x(8)),-cos(x(7))*sin(x(8))];

%          x y z vx     vy    vz  r p  t
df_perp_mat      = @(x) [zeros(3), eye(3),zeros(3,3);
    zeros(3,6), -(db_T_q(x)*[1;0])*x(9), -(db_T_q(x)*[0;1])*x(9), -b_T(x)];

eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N, ctrl_N,2);
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);

for i4 = 1:length(r_range)
    for i5 = 1:length(p_range)
        for i6 = 1:length(th_range)
            
            x = [zeros(6,1);
                r_range(i4); p_range(i5); th_range(i6)];
            
            W = W_eval(w_poly_fnc(x));
            M = W\eye(n);
            Theta = chol(M);
            Theta_Bw = Theta*Bw;
            sigma_ThBw(i4,i5,i6) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));
            
            L = chol(W);
            
            f = f_mat(x);
            df_perp = df_perp_mat(x);
            F = df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W(1:6,1:6);
            
            R_CCM = -F;
            
            eig_CCM(i4,i5,i6) = min(eig(R_CCM));
            eig_W(i4,i5,i6,1) = min(eig(W));
            eig_W(i4,i5,i6,2) = max(eig(W));
        end
    end
end

d_bar = max(sigma_ThBw(:))/lambda;
disp('d_bar'); disp(d_bar);
disp('W:'); disp(min(min(min(eig_W(:,:,:,1)))));
disp(max(max(max(eig_W(:,:,:,2)))));
disp('CCM:'); disp(min(eig_CCM(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper_mat)));






