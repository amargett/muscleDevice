clear; close all; 
global data; global data2; global p; global x_vals; global y_vals; 
data = []; data2 = []; p = []; x_vals = []; y_vals = []; 

%% Parameters to set
E = 1.685E9; % material property
r_well = 7.5E-3; % mm
t0 = 0.25E-3; % mm, thickness of constant part of beam
r_inner = 1E-3; % mm, radius of inner circle
l_m = 1E-3; % minimum space for pushing motor
SF = 5; % motor torque safety factor
Tmax = 0.2292; % Nm
sigma_yield = 0.5E9; % Pa, material yield stress
poly_order = 5;
theta_f = 7*pi/8; 

params = [E, r_well, t0, r_inner, poly_order, SF, Tmax, sigma_yield, theta_f];

%% Constraints
k_min = 0.3; k_max = 5; %N/m
l_max = r_well- t0 - r_inner - l_m; % upper bound beam length
% manufacturing limits, not sure of this yet
t_min = 0.1E-3; %mm, min reasonable thickness
w_min = 0.1E-3; %mm, min reasonable width

%% Optimization to minimize Kb/Kr
fun = @(z) get_prbm_cost(z, params); 

lb = [3*pi/8, r_well/2, t_min, w_min]; % lower bound 
ub = [80*pi/180, l_max, 0.5E-3, 0.5E-3]; % upper bound
z0 = [0, r_well/2, t_min, w_min]; 

[z, fval1] = fmincon(fun, z0, [], [], [], [], lb, ub); 

z_deg_mm = z.*[180/pi 1E3 1E3 1E3]; 

get_prbm_cost(z, params); %update data one more time

%% Optimization to minimize Klow/Khigh
p_beams = [E t0 r_well z(4)]; 
fun = @(h) get_beam_cost(h, p_beams, p); 
h0= [-0.006879 0.07310649 -0.29553705 0.40500924 -0.32101599 5.89045809];
% 
ms = MultiStart('UseParallel', true); % Enable parallel computing for speed
problem = createOptimProblem('fmincon', 'x0', h0, 'objective', fun, 'lb', [], 'ub', []);
[x_multi, fval_multi] = run(ms, problem, 50); % Try 50 different initial points

h = x_multi; 
fval2 = fval_multi;
fun(h)
% h = h0; 
% fval2 = fun(h);
%% Display Variables & Data
varNames = {'theta [deg]', 'l [mm]', 't [mm]', 'w [mm]'};
resultsTable = table(varNames', z_deg_mm', 'VariableNames', ...
    {'Variable', 'Optimal_Value'});
disp(resultsTable);
disp('Minimum Ratio Kb/Kr: '); disp(fval1); 


Kb = data(1); Kr = data(2); Kb_motor = data(3); Kb_yield = data(4); 
K_high = data2(1); K_low = data2(2); 

disp('Bending Stiffness Kb: '); disp(Kb);
disp('Radial Stiffness Kr: '); disp(Kr); 
disp('Kb Motor: '); disp(Kb_motor);
disp('Kb Yield: '); disp(Kb_yield);

disp('Beam Stiffness Ratio: '); disp(fval2);
disp('K High: '); disp(K_high); 
disp('K Low: '); disp(K_low); 

%% Create & Display PRBM Polynomial String 
poly_eq = sprintf('%.6f*x^%d', p(1), poly_order); % Start with the highest-order term
for k = 2:length(p)
    poly_eq = sprintf('%s + %.8f*x^%d', poly_eq, p(k), poly_order - (k - 1));
end

% Evaluate the polynomial for a smooth curve
x_fit = linspace(min(x_vals), max(x_vals), 100);
y_fit = polyval(p, x_fit);

y_fit2 = polyval(h, x_fit); 

% Plot original data
figure;
plot(x_vals, y_vals, 'o', 'DisplayName', 'Original Data');
hold on;

% Plot polynomial fit
plot(x_fit, y_fit, '-', 'DisplayName', ['PRBM Poly Fit (Order ', num2str(poly_order), ')']);

% Plot polynomial fit
plot(x_fit, y_fit2, '-', 'DisplayName', ['Beam Poly Fit (Order ', num2str(poly_order), ')']);
hold off;

% Add labels and legend
xlabel('X Values');
ylabel('Y Values');
title('Data and Polynomial Fit');
legend('Location', 'Best');
grid on;

disp('PRBM polynomial: '); disp(poly_eq)

%% Create & Display Beam Polynomial String 
poly_eq = sprintf('%.6f*x^%d', h(1), poly_order); % Start with the highest-order term
for k = 2:length(h)
    poly_eq = sprintf('%s + %.8f*x^%d', poly_eq, h(k), poly_order - (k - 1));
end
disp('Beam polynomial: '); disp(poly_eq)

