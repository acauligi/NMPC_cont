function dx_dot = quad_sim_geom(t,x,...
    t_vec,state_nom,ctrl_nom,ang_a,f,B,B_w,w)

x_nom = interp1(t_vec,state_nom,t); x_nom = x_nom';
u_nom = interp1(t_vec,ctrl_nom,t); u_nom = u_nom';
om_d_nom = interp1(t_vec,ang_a,t); om_d_nom = om_d_nom';

m = 4.34;
g = 9.81;
Jx = 0.082; Jy =  0.0845; Jz =  0.1377;
Jq = diag([Jx;Jy;Jz]); 

kx = 16*m;
kv = 5.6*m;
kR = 8.81;
k_om = 2.54;

B_des = B(x_nom);
a_des = B_des(4:6,1)*u_nom(1) + [0;0;g];

R = rot_matrix(x(9),x(8),x(7));
R_des = rot_matrix(x_nom(9),x_nom(8),x_nom(7));

ex = x(1:3) - x_nom(1:3);
ev = x(4:6) - x_nom(4:6);

thrust = (-kx*ex - kv*ev - m*g*[0;0;1] + m*a_des)'*(R*[0;0;-1]);

eR = vee(0.5*(R_des'*R - R'*R_des));

om = x(10:12);
om_des = x_nom(10:12);

e_om = om - R'*R_des*om_des;

torque = -kR*eR - k_om*e_om + Skew(om)*Jq*om -...
         Jq*(Skew(om)*R'*R_des*om_des - R'*R_des*om_d_nom);
     
u = [thrust; torque];     

dx_dot = f(x) + B(x)*u + B_w(x)*w;

end