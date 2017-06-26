load_PVTOL_config;
load MP_WARM_PVTOL.mat;

ts = 0.5*mp_warm.Tp*(mp_warm.s_t+1);

A = [zeros(3) eye(3); zeros(3,6)];
B = [zeros(3); eye(3)];
C = eye(6);
D = zeros(6,3);
system = ss(A,B,C,D);

ts = linspace(0,mp_warm.Tp,numel(mp_warm.s_t));
output = lsim(system, mp_warm.ctrl', ts, test_state);
mp_warm.state = output;
save('MP_WARM_PVTOL.mat','mp_warm');