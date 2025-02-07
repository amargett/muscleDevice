global data; global p; 
data = []; p = []; params = []; 

%Parameters to set
E = 18E9; % material property
r_well = 7.5E-3; 
t0 = 1E-3; % thickness of constant part of beam
r_inner = 0.5E-3; % radius of inner circle
l_tip = 0.5E-3; % length of variable stiffness tip
l_m = 1E-3; % minimum space for pushing motor
SF = 5; %motor torque safety factor
Tmax = 0.2292; % Nm
sigma_yield = 5E7; % Pa, material yield stress

params = [E, r_well, t0, r_inner, l_tip, SF, Tmax, sigma_yield];

% Constraints
k_min = 0.3; k_max = 5; %N/m
l_max = r_well- t0 - r_inner - l_m - l_tip; % upper bound beam length
% manufacturing limits, not sure of this yet
t_min = 0.1E-3; %m, min reasonable thickness
w_min = 0.1E-3; %m, min reasonable width


% optimization to minimize Kb/Kr
fun = @(z) get_prbm_cost(z, params); 

lb = [0, r_well/2, t_min, w_min]; % lower bound 
ub = [pi/2, l_max, 0.5E-3, 0.5E-3]; % upper bound
z0 = [0, r_well/2, t_min, w_min]; 

[z, fval] = fmincon(fun, z0, [], [], [], [], lb, ub); 
z_deg_mm = z.*[180/pi 1E3 1E3 1E3]; 

varNames = {'theta [deg]', 'l [mm]', 't [mm]', 'w [mm]'};
resultsTable = table(varNames', z_deg_mm', 'VariableNames', ...
    {'Variable', 'Optimal_Value'});
disp(resultsTable);
disp('Minimum Ratio Kb/Kr: '); disp(fval); 

get_prbm_cost(z, params); %update data one more time
Kb = data(1); Kr = data(2); Kb_motor = data(3); Kb_yield = data(4); s = data(5); 
disp('Bending Stiffness Kb: '); disp(Kb);
disp('Radial Stiffness Kr: '); disp(Kr); 
disp('Kb Motor: '); disp(Kb_motor);
disp('Kb Yield: '); disp(Kb_yield);
disp('S: '); disp(s);


poly_order = 4; 
poly_eq = sprintf('y = %.6fx^%d', p(1), poly_order); % Start with the highest-order term
for k = 2:length(p)
    poly_eq = sprintf('%s + %.8f*x^%d', poly_eq, p(k), poly_order - (k - 1));
end

disp('PRBM polynomial: '); disp(poly_eq)



