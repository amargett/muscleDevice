clear; close all; 
global data; global data2; global p; global x_vals; global y_vals; 
data = []; data2 = []; p = []; x_vals = []; y_vals = []; 

%% Parameters to set
E = 0.14E9; % material property
r_well = 7.5E-3; % mm
r_inner = 1E-3; % mm, radius of inner circle
l_m = 1E-3; % minimum space for pushing motor
SF = 5; % motor torque safety factor
Tmax = 0.2292; % Nm
sigma_yield = 0.5E9; % Pa, material yield stress
poly_order = 5;
theta_f = 7*pi/8; 
w = 3.175E-3; 
%% Find K_low, K_high numbers
t0_vals = 10^-3 *[0.25 0.5 0.75 1.0];

for t0 = t0_vals
    k_min = 0.3; k_max = 5; %N/m
    l_max = r_well- t0 - r_inner - l_m; % upper bound beam length
    % manufacturing limits, not sure of this yet
    t_min = 1E-3; %mm, min reasonable thickness
    w_min = 3.175E-3; %mm, min reasonable width

    params = [E, r_well, t0, r_inner, poly_order, SF, Tmax, sigma_yield, theta_f];

    fun = @(z) get_prbm_cost(z, params); 
    
    lb = [3*pi/8, r_well/2, 0.1E-3, w_min]; % lower bound 
    ub = [80*pi/180, l_max, 0.5E-3, 3.175E-3]; % upper bound
    z0 = [0, r_well/2, t_min, w_min]; 
    
    [z, fval1] = fmincon(fun, z0, [], [], [], [], lb, ub); 
    
    z_deg_mm = z.*[180/pi 1E3 1E3 1E3]; 
    
    get_prbm_cost(z, params); %update data one more time
    
    poly = p;
    
    N= 10; 
    % t_min_vals = linspace(t_min, t0, N); 
    % 
    % K_high_vals = zeros(size(t_min_vals)); K_low_vals = zeros(size(t_min_vals));
    % 
    % for i = 1:N
    %     [K_high_vals(i), K_low_vals(i)] = get_beam_k(poly, E, t0, r_well, w, t_min_vals(i)); 
    % end

    theta_b_vals = linspace(pi/2, pi, N); 
    K_vals = zeros(size(theta_b_vals));

    for i = 1:N
        K_vals(i) = get_beam_k(poly, E, t0, r_well, w, t_min, theta_b_vals(i)); 
    end
    plot(theta_b_vals*180/pi, K_vals); 
    hold on

end

xlabel('theta b [deg]'); ylabel('K beam [N/m]'); 
legend('t0 = 0.25 mm', 't0 = 0.5 mm', 't0 = 0.75 mm', 't0 = 1 mm');

hold off