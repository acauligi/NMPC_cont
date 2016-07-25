clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
%% Constants

n = 14;
m = 4;

Jx = 0.082; Jy =  0.0845; Jz =  0.1377;
Jq = diag([Jx;Jy;Jz]);
a = (Jy-Jz)/Jx; b = (Jz-Jx)/Jy; c = (Jx-Jy)/Jz;
mq = 4.34;
g = 9.81;

%% Generate desired trajectory

T = 10;
dt = 0.001;
t = (0:dt:T)';

[state_nom, ctrl_nom, ang_a] = generate_quad_traj(t,Jq,mq,g);

%Plotting
figure()
subplot(3,1,1)
plot(t,state_nom(:,9:11)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Attitude [deg]');
legend('Roll','Pitch','Yaw');

subplot(3,1,2)
plot(t,state_nom(:,12:14)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Angular rate [deg/s]');
legend('p','q','r');

subplot(3,1,3)
plot(t,ang_a(:,1:3)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Angular rate [deg/s^2]');
legend('pd','qd','rd');

figure()
subplot(3,1,1)
plot(t,state_nom(:,7),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('T [N]');

subplot(3,1,2)
plot(t,state_nom(:,8),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('$\dot{T}$ [N/s]','interpreter','latex');

subplot(3,1,3)
plot(t,ctrl_nom(:,1),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('$\ddot{T}$ [N/s]','interpreter','latex');

%% Construct constant CCM

lambda = 1.0;
r_vec = [4;4;4;2];
Ac = [];
Bc = [];
for i = 1:m
    Ac_i = [zeros(r_vec(i)-1,1), eye(r_vec(i)-1);
            zeros(1,r_vec(i))];
    Ac = blkdiag(Ac,Ac_i);
    
    Bc_i = [zeros(r_vec(i)-1,1);1];
    Bc = blkdiag(Bc,Bc_i);
end
B_perp = null(Bc');
cvx_begin sdp
    variable W_ccm(n,n) symmetric
    variables w_lower w_upper
    minimize (w_upper - w_lower)
    subject to
    W_ccm >= w_lower*eye(n);
    W_ccm <= w_upper*eye(n);
    w_lower >= 0.01;

    B_perp'*(Ac*W_ccm + W_ccm*Ac' + 2*lambda*W_ccm)*B_perp <= 0;
cvx_end

M_ccm = W_ccm\eye(n);

pause;

%% Setup aux controller and dynamics

eps_u = 0.1;
aux_prob = setup_opt_aux(m,eps_u);
geo_Ke = 1;

% dq = 0.2;
%x:[x,y,z,vx,vy,vz,th,th_d,phi,th,psi,p,q,r]

roll_d = @(x) x(12) + x(13)*tan(x(10))*sin(x(9)) + x(14)*tan(x(10))*cos(x(9));
pitch_d = @(x) x(13)*cos(x(9)) - x(14)*sin(x(9));
yaw_d = @(x) x(13)*sec(x(10))*sin(x(9)) + x(14)*sec(x(10))*cos(x(9));

bq_1 = @(x) sin(x(9))*sin(x(11)) + cos(x(9))*sin(x(10))*cos(x(11));
bq_2 = @(x) -sin(x(9))*cos(x(11)) + cos(x(9))*sin(x(10))*sin(x(11));
bq_3 = @(x) cos(x(9))*cos(x(10));

f = @(x) [x(4); 
          x(5); 
          x(6);
          -(1/mq)*x(7)*bq_1(x);
          -(1/mq)*x(7)*bq_2(x);
          g - (1/mq)*x(7)*bq_3(x);
          x(8);
          0;
          roll_d(x);
          pitch_d(x);
          yaw_d(x);
          a*x(13)*x(14);
          b*x(12)*x(14);
          c*x(12)*x(13)];
          
B = [zeros(7,4);
     1, zeros(1,3);
     zeros(3,4);
     zeros(3,1), (Jq\eye(3))];

B_w = B;
w_max = max(norms(0.1*ctrl_nom(:,2:4),2,2));
sigma_Bw = sqrt(max(eig(B_w'*B_w)));
d_bar = sigma_Bw*w_max*sqrt(max(eig(M_ccm)))/lambda;
% euc_bound = d_bar/sqrt(min(eig(M_ccm)));

%% Setup geodesic mapping

phi = @(x) [x(1);
            x(4);
            -(1/mq)*x(7)*bq_1(x);
            -(1/mq)*x(8)*bq_1(x)-(1/mq)*x(7)*(roll_d(x)*cos(x(9))*sin(x(11))+yaw_d(x)*sin(x(9))*cos(x(11))-roll_d(x)*sin(x(9))*sin(x(10))*cos(x(11))+pitch_d(x)*cos(x(9))*cos(x(10))*cos(x(11))-yaw_d(x)*cos(x(9))*sin(x(10))*sin(x(11)));
            x(2);
            x(5);
            -(1/mq)*x(7)*bq_2(x);
            -(1/mq)*x(8)*bq_2(x)-(1/mq)*x(7)*(-roll_d(x)*cos(x(9))*cos(x(11))+yaw_d(x)*sin(x(9))*sin(x(11))-roll_d(x)*sin(x(9))*sin(x(10))*sin(x(11))+pitch_d(x)*cos(x(9))*cos(x(10))*sin(x(11))+yaw_d(x)*cos(x(9))*sin(x(10))*cos(x(11)));
            x(3);
            x(6);
            g-(1/mq)*x(7)*bq_3(x);
             -(1/mq)*x(8)*bq_3(x)-(1/mq)*x(7)*(-roll_d(x)*sin(x(9))*cos(x(10))-pitch_d(x)*cos(x(9))*sin(x(10)));
            x(11);
            yaw_d(x)];

phi_d = @(x) [ 1, 0, 0, 0, 0, 0,                                                                                  0,                                                  0,                                                                                                                    0,                                                                               0,                                                                                                                                                0,                                                       0,                          0,                0;
               0, 0, 0, 1, 0, 0,                                                                                  0,                                                  0,                                                                                                                    0,                                                                               0,                                                                                                                                                0,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 0,                                 -(sin(x(9))*sin(x(11)) + cos(x(9))*cos(x(11))*sin(x(10)))/mq,                                                  0,                                                              -(x(7)*(cos(x(9))*sin(x(11)) - cos(x(11))*sin(x(9))*sin(x(10))))/mq,                                              -(x(7)*cos(x(9))*cos(x(10))*cos(x(11)))/mq,                                                                                          -(x(7)*(cos(x(11))*sin(x(9)) - cos(x(9))*sin(x(10))*sin(x(11))))/mq,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 0, -(x(13)*cos(x(10))*cos(x(11)) + x(12)*cos(x(9))*sin(x(11)) - x(12)*cos(x(11))*sin(x(9))*sin(x(10)))/mq, -(sin(x(9))*sin(x(11)) + cos(x(9))*cos(x(11))*sin(x(10)))/mq, (x(7)*x(12)*sin(x(9))*sin(x(11)) - x(8)*cos(x(9))*sin(x(11)) + x(8)*cos(x(11))*sin(x(9))*sin(x(10)) + x(7)*x(12)*cos(x(9))*cos(x(11))*sin(x(10)))/mq, (cos(x(11))*(x(7)*x(13)*sin(x(10)) - x(8)*cos(x(9))*cos(x(10)) + x(7)*x(12)*cos(x(10))*sin(x(9))))/mq, -(x(8)*cos(x(11))*sin(x(9)) + x(7)*x(12)*cos(x(9))*cos(x(11)) - x(7)*x(13)*cos(x(10))*sin(x(11)) - x(8)*cos(x(9))*sin(x(10))*sin(x(11)) + x(7)*x(12)*sin(x(9))*sin(x(10))*sin(x(11)))/mq, -(x(7)*(cos(x(9))*sin(x(11)) - cos(x(11))*sin(x(9))*sin(x(10))))/mq, -(x(7)*cos(x(10))*cos(x(11)))/mq,                0;
               0, 1, 0, 0, 0, 0,                                                                                  0,                                                  0,                                                                                                                    0,                                                                               0,                                                                                                                                                0,                                                       0,                          0,                0;
               0, 0, 0, 0, 1, 0,                                                                                  0,                                                  0,                                                                                                                    0,                                                                               0,                                                                                                                                                0,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 0,                                  (cos(x(11))*sin(x(9)) - cos(x(9))*sin(x(10))*sin(x(11)))/mq,                                                  0,                                                               (x(7)*(cos(x(9))*cos(x(11)) + sin(x(9))*sin(x(10))*sin(x(11))))/mq,                                              -(x(7)*cos(x(9))*cos(x(10))*sin(x(11)))/mq,                                                                                          -(x(7)*(sin(x(9))*sin(x(11)) + cos(x(9))*cos(x(11))*sin(x(10))))/mq,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 0,  (x(12)*cos(x(9))*cos(x(11)) - x(13)*cos(x(10))*sin(x(11)) + x(12)*sin(x(9))*sin(x(10))*sin(x(11)))/mq,  (cos(x(11))*sin(x(9)) - cos(x(9))*sin(x(10))*sin(x(11)))/mq, (x(8)*cos(x(9))*cos(x(11)) - x(7)*x(12)*cos(x(11))*sin(x(9)) + x(8)*sin(x(9))*sin(x(10))*sin(x(11)) + x(7)*x(12)*cos(x(9))*sin(x(10))*sin(x(11)))/mq, (sin(x(11))*(x(7)*x(13)*sin(x(10)) - x(8)*cos(x(9))*cos(x(10)) + x(7)*x(12)*cos(x(10))*sin(x(9))))/mq, -(x(8)*sin(x(9))*sin(x(11)) + x(7)*x(13)*cos(x(10))*cos(x(11)) + x(7)*x(12)*cos(x(9))*sin(x(11)) + x(8)*cos(x(9))*cos(x(11))*sin(x(10)) - x(7)*x(12)*cos(x(11))*sin(x(9))*sin(x(10)))/mq,  (x(7)*(cos(x(9))*cos(x(11)) + sin(x(9))*sin(x(10))*sin(x(11))))/mq, -(x(7)*cos(x(10))*sin(x(11)))/mq,                0;
               0, 0, 1, 0, 0, 0,                                                                                  0,                                                  0,                                                                                                                    0,                                                                               0,                                                                                                                                                0,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 1,                                                                                  0,                                                  0,                                                                                                                    0,                                                                               0,                                                                                                                                                0,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 0,                                                             -(cos(x(9))*cos(x(10)))/mq,                                                  0,                                                                                             (x(7)*cos(x(10))*sin(x(9)))/mq,                                                        (x(7)*cos(x(9))*sin(x(10)))/mq,                                                                                                                                                0,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 0,                                           (x(13)*sin(x(10)) + x(12)*cos(x(10))*sin(x(9)))/mq,                             -(cos(x(9))*cos(x(10)))/mq,                                                                          (cos(x(10))*(x(8)*sin(x(9)) + x(7)*x(12)*cos(x(9))))/mq,            (x(7)*x(13)*cos(x(10)) + x(8)*cos(x(9))*sin(x(10)) - x(7)*x(12)*sin(x(9))*sin(x(10)))/mq,                                                                                                                                                0,                                (x(7)*cos(x(10))*sin(x(9)))/mq,           (x(7)*sin(x(10)))/mq,                0;
               0, 0, 0, 0, 0, 0,                                                                                  0,                                                  0,                                                                                                                    0,                                                                               0,                                                                                                                                                1,                                                       0,                          0,                0;
               0, 0, 0, 0, 0, 0,                                                                                  0,                                                  0,                                                                                 (x(13)*cos(x(9)) - x(14)*sin(x(9)))/cos(x(10)),             (sin(x(10))*(x(14)*(2*sin(x(9)/2)^2 - 1) - x(13)*sin(x(9))))/(sin(x(10))^2 - 1),                                                                                                                                                0,                                                       0,           sin(x(9))/cos(x(10)), cos(x(9))/cos(x(10))];

M = @(x) M_ccm*phi_d(x);

%% Continuous Simulation

ode_options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
start_p = state_nom(1,:);

% [t_vec,x_act] = ode45(@(t_vec,x_act)quad_sim_cont(t_vec,x_act,...
%         t,state_nom,ctrl_nom,phi,M_ccm,M,lambda,f,B,B_w,w_dist),...
%         t,start_p,ode_options);

%% Discrete Simulation

T_steps = length(t)-1;

t_opt = cell(T_steps,1);

x_act = zeros(T_steps+1,n);
x_act(1,:) = start_p;

state = cell(T_steps,1);

ctrl = zeros(T_steps,m);
aux_ctrl = zeros(T_steps,m);

solved = ones(T_steps,1);

E = zeros(T_steps,1);

u_prev = zeros(m,1);

for i = 1:T_steps
    
%     fprintf('%d/%d \n',i, T_steps);
    
    x_nom = state_nom(i,:)';
    u_nom = ctrl_nom(i,:)';
    
    xi_nom = phi(x_nom);
    xi_act = phi(x_act(i,:)');
    
    X_dot = kron(ones(1,2),xi_act-xi_nom);
    X = [x_nom, x_act(i,:)'];
    E(i) = (xi_act - xi_nom)'*M_ccm*(xi_act - xi_nom);
    
    [aux, solved(i)] = compute_opt_aux_FL(aux_prob,geo_Ke,...
                            X,X_dot,E(i),M,f,B,u_nom,u_prev,eps_u,lambda);
    aux_ctrl(i,:) = aux';
    
    ctrl(i,:) = u_nom' + aux';%zeros(1,m);
    u_prev = ctrl(i,:)';
    
    w_dist = u_nom.*[0;0.1*ones(3,1)];
    x_act(i+1,:) = x_act(i,:)+(f(x_act(i,:)') + B*ctrl(i,:)' + B_w*w_dist)'*dt;
    
end


%% Plot

close all

%Trajectory plot
figure()
plot3(x_act(:,1),x_act(:,2),x_act(:,3),'b-','linewidth',2); hold on
plot3(state_nom(:,1),state_nom(:,2),state_nom(:,3),'r-','linewidth',2);
grid on
xlabel('x'); ylabel('y'); zlabel('h');
set(gca,'ZDir','Reverse');
set(gca,'YDir','Reverse');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Trajectory errors
figure()
plot(t,state_nom(:,1:3)-x_act(:,1:3),'linewidth',2);
grid on
xlabel('Time [s]');
ylabel('Traj errors');
legend('e_x','e_y','e_z');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%Control effort
figure()
subplot(3,1,1)
plot(t(1:end-1),ctrl_nom(1:end-1,1),'r-','linewidth',2); hold on
plot(t(1:end-1),ctrl(:,1),'b-','linewidth',2);
xlabel('Time [s]');
ylabel('$\ddot{u}_1$','interpreter','latex'); 
grid on
legend('nominal','net');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,2)
plot(t,state_nom(:,7),'r-','linewidth',2); hold on
plot(t,x_act(:,7),'b-','linewidth',2);
xlabel('Time [s]');
ylabel('Thrust [N]'); 
grid on
legend('nominal','net');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,3)
plot(t(1:end-1),ctrl_nom(1:end-1,2:4),'-','linewidth',2); hold on
plot(t(1:end-1),ctrl(:,2:4),'--','linewidth',2);
xlabel('Time [s]'); 
ylabel('Torque [Nm]'); 
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


% Geodesic Energy
figure()
plot(t(1:end-1),E,'b-','linewidth',2); hold on
plot(t(1:end-1),(d_bar^2)*ones(T_steps,1),'r-','linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Energy');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

% Solve success
figure()
plot(t(1:end-1),solved(:,1),'go','markersize',10,'markerfacecolor','g');
grid on
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)



%% Animate

plot_quad_movie(state_nom(:,1),state_nom(:,2),state_nom(:,3),t,x_act,20,n)

