clear all; close all; clc;
yalmip('clear');

%% Allgower example

lambda_range = 1.742857142857143; %final optimal
% lambda_range = linspace(0.8,2,15); %line search range
bounds = nan(15,1);
eps = 1e-3;
for ll = 1:length(lambda_range) 
    lambda = lambda_range(ll);
    fprintf('lambda: %f ', lambda);
    cond_l = 1.6336;
    cond_u = 1.64;
    cond_best = cond_u;
    
    while (cond_u-cond_l>eps) %bisection search on condition number
        condn = (cond_l+cond_u)/2;
        
        fprintf(' cond: %f ', condn);
        
        yalmip('clear');
        w_upper = sdpvar(1); w_lower = sdpvar(1);
        
        x1 = sdpvar(1);
        x2 = sdpvar(1);
        
        x = [x1 x2]';
        
        f = [-1*x1 + 2*x2;
            -3*x1 + 4*x2 - 0.25*(x2^3)];
        B = [0.5;-2];
        df = jacobian(f,x);
               
        delta = sdpvar(2,1);
        
        [W_11,c_11]= polynomial(x,4,0);
        [W_12,c_12]= polynomial(x,4,0);
        [W_22,c_22]= polynomial(x,4,0);
        
        [rho,c_rho] = polynomial(x,4,0);
        
        W = [W_11,W_12;
            W_12,W_22];
        
        dW_f = jacobian(W(:),x)*f;
        dW_f = reshape(dW_f,2,2);
        
        %Box constraints
        g = [x1+5; 5-x1; x2+5; 5-x2];
        
        %Box Lagrange functions
        
        %Upper and Lower definiteness
        [su_1,c_us_1] = polynomial([x;delta],4,0);
        [su_2,c_us_2] = polynomial([x;delta],4,0);
        [su_3,c_us_3] = polynomial([x;delta],4,0);
        [su_4,c_us_4] = polynomial([x;delta],4,0);
        [sl_1,c_ls_1] = polynomial([x;delta],4,0);
        [sl_2,c_ls_2] = polynomial([x;delta],4,0);
        [sl_3,c_ls_3] = polynomial([x;delta],4,0);
        [sl_4,c_ls_4] = polynomial([x;delta],4,0);
        
        %CCM condition
        [s_ccm_1,c_ccm_1] = polynomial([x;delta],4,0);
        [s_ccm_2,c_ccm_2] = polynomial([x;delta],4,0);
        [s_ccm_3,c_ccm_3] = polynomial([x;delta],4,0);
        [s_ccm_4,c_ccm_4] = polynomial([x;delta],4,0);
        
        %CCM condition
        R_CCM = -(-dW_f + df*W + W*df' - rho*(B*B') + 2*lambda*W);
        p_CCM = delta'*R_CCM*delta - [s_ccm_1,s_ccm_2,s_ccm_3,s_ccm_4]*g;        
        
        %Killing field condition
        dW_B = jacobian(W(:),x)*B;
        dW_B = reshape(dW_B,2,2);
        killing_11 = dW_B(1,1);
        killing_11_2 = (-dW_B(1,1));
        
        killing_12 = dW_B(1,2);
        killing_12_2 = (-dW_B(1,2));
        
        killing_22 = dW_B(2,2);
        killing_22_2 = (-dW_B(2,2));
        
        %W uniform bounds
        W_bounds = [w_lower>=0.0035; w_upper>=w_lower];
        
        %Condition bound
        W_cond = [w_upper-condn*w_lower <= 0];
        
        %W pos def
        p_W_up = (w_upper*(delta'*delta)-delta'*W*delta) - ...
            [su_1,su_2,su_3,su_4]*g;
        
        p_W_low = (delta'*W*delta - w_lower*(delta'*delta))-...
            [sl_1,sl_2,sl_3,sl_4]*g;
        
        coeff_List = [c_rho;c_11;c_12;c_22;
            c_us_1;c_us_2;c_us_3;c_us_4;c_ls_1;c_ls_2;c_ls_3;c_ls_4;
            c_ccm_1;c_ccm_2;c_ccm_3;c_ccm_4];
        
        constr_List = [ sos(p_W_low);sos(p_W_up);
            sos(p_CCM);
            sos(su_1); sos(su_2);
            sos(su_3); sos(su_4);
            sos(sl_1); sos(sl_2);
            sos(sl_3); sos(sl_4);
            W_bounds;W_cond;
            sos(killing_11); sos(killing_11_2);
            sos(killing_12); sos(killing_12_2);
            sos(killing_22); sos(killing_22_2);
            sos(s_ccm_1); sos(s_ccm_2); sos(s_ccm_3); sos(s_ccm_4)];
        
        options = sdpsettings('solver','mosek','verbose',0);
        SOS_soln = solvesos(constr_List,[],options,coeff_List);
        if (SOS_soln.problem == 0)
            fprintf('feasible\n');
            cond_u = condn;
            cond_best = condn;
            coeff_sol = clean(double(coeff_List),1e-3);
            W_sol = replace(W,coeff_List,coeff_sol);
            dW_f = jacobian(W(:),x)*f; dW_f = reshape(dW_f,2,2);
            dW_f_sol = replace(dW_f,coeff_List,coeff_sol);
            rho_sol = replace(rho,coeff_List,coeff_sol);
        else
            cond_l = condn;
        end
    end
    fprintf('\n');
    bounds(ll) = (sqrt(condn))/lambda;
end
pause;
clc;

% figure()
% plot(lambda_range, bounds,'ro','markerfacecolor','g','markersize',20);
% grid on
% xlabel('\lambda');
% ylabel('$\|x^{*}-x\|/\bar{w}$','Interpreter','Latex');
% % title('Robustness optimization');
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
% % Display
% 
% disp('W')
% sdisplay(W_sol)
% 
% disp('rho')
% sdisplay(rho_sol)
% 
% disp('Bound achieved');
% disp((1/lambda)*sqrt(condn));

%% Compute control bounds using optimal metric (W = const matrix)

M = W_sol\eye(2);
alpha = max(eig(M));
w = 0.1;
d_bar = sqrt(alpha)*w/lambda;

%Compare RCI sets
P_rci = diag([39.0251, 486.0402]);

figure()
Ellipse_plot(M*(1/d_bar^2),[0;0],20,'k');
Ellipse_plot(P_rci,[0;0],20,'r');

L = chol(W_sol);
F = df*W_sol + W_sol*df' + 2*lambda*W_sol;

x2_range = linspace(-5,5,30);
delta_u = zeros(length(x2_range),1);

for i = 1:length(x2_range)
    F_sol = replace(F,x2,x2_range(i));
    
    delta_u(i) = 0.5*d_bar*max(eig((inv(L))'*F_sol*inv(L)))/...
                          sqrt(max(eig((inv(L))'*(B*B')*inv(L))));
end

figure()
plot(x2_range,delta_u,'b-','linewidth',2);
grid on
xlabel('x2'); ylabel('$\bar{\delta}_u$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)


