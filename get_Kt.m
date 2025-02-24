function Kt = get_Kt(params)
% params: [r_inner, r_well, t0, l_m, t_min, w, E]
r_inner = params(1); r_well = params(2); t0 = params(3); l_m = params(4); 
t_min = params(5); w = params(6); l_tip = params(7); E = params(8); 

l = r_well -t0-r_inner -l_m -l_tip; 
t = t_min; 
I = @(t,w) w*t^3/12; 
gamma = 0.8517; K_theta = 2.65; % PRBM Coefficients

% Calculating stiffness for two beams
Kt = 4*gamma*K_theta*E*I(t,w)/l % PRBM eqn, dT/dtheta

delta_FEA= 3.615E-5; F_FEA = 3E-7; 

K_FEA = F_FEA*(l+l_m+l_tip)*r_well/delta_FEA; 
disp('K found in FEA: '); disp(K_FEA);
end