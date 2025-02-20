function cost = get_prbm_cost(z, params)
    % Initiliazing global variables to track Kb, Kr, p
    global data
    global p
    global x_vals
    global y_vals
    
    %% Get Input Variables
    E = params(1); r_well = params(2); t0 = params(3); r_inner = params(4); 
    poly_order = params(5); SF = params(6); Tmax = params(7); sigma_yield = params(8); 
    theta_f = params(9);

    theta = z(1); l = z(2); t = z(3); w = z(4); 
    
    %% PRBM/geometric parameters
    gamma = 0.8517; 
    l_r = r_well-r_inner-t0-l-0.75E-3; % rigid length
    l_tip = r_well-t0-(l+l_r+r_inner)*sin(theta); % mm, length of variable stiffness tip
    lc=(1-gamma)*l/2;
    l_top = (r_inner+l+l_r)*cos(theta); 

    %% Fit Polynomial to PRBM Data
    N = 100; theta_vals = linspace(0, theta, N); 
    
    rtipdata = zeros(2,N); 
    for i = 1:N
        rtipdata(:,i) = prbm(theta_vals(i)); 
    end
    
        % Given x value to search for
    x_target = 0.844E-3; % Change this to x value from FEA
    
    % Find the index of the closest x value in the first row of rtipdata
    [~, idx] = min(abs(rtipdata(1, :) - x_target));
    
    % Retrieve the corresponding theta value
    theta_closest = theta- theta_vals(idx);

    fea_K = (l+l_r-lc)*0.1/theta_closest; 
   
    % Display result
    disp('K Value found in FEA: '+ string(fea_K)); 
    
    % Extract x and y components
    x_vals = rtipdata(1, :)*10^3; y_vals = rtipdata(2, :)*10^3;

    % Fit a polynomial to the data
    p = polyfit(x_vals, y_vals, poly_order);
    
    %% Radial Stiffness Kr
    I = @(t,w) w*t^3/12; 
    k1 = (l*(sin(theta))^2/(E*w*t));
    k2 = (cos(theta))^2*l^3/(3*E*I(t,w));
    Kr = 1/(l*(sin(theta))^2/(E*w*t)+ (cos(theta))^2*l^3/(3*E*I(t,w)));
    
    %% Bending Stiffness Kb
    K_theta = 2.65; % PRBM coefficient
    Kb = 4*gamma*K_theta*E*I(t,w)/l; % PRBM eqn, dT/dtheta

    % Kb Upper Bound due to Motor Specs
    Kb_max_motor = 2*Tmax/(SF*pi/2); 
    
    % Kb Upper Bound due to yield stress;
    Kb_max_yield = 2*sigma_yield*w*t^2/(6*pi/4); 
    
    % penalizes heavily for breaking bounds
    penalty = 100 * (Kb > Kb_max_motor | Kb > Kb_max_yield);
    
    cost = Kb/Kr + penalty; 
    data = [Kb Kr Kb_max_motor Kb_max_yield];

    %% Helper Functions
    function rC = prbm(theta_p)
        % finds endpoint position as a function of theta
        ihat = [1 0]; jhat = [0 1]; 
        theta2 = theta_p + pi - theta; 
        theta3 = theta_p + pi/2 -theta; 
        rA = (lc+ gamma*l*sin(theta_p))*jhat + gamma*l*cos(theta_p)*ihat; 
        rC = rA + l_top/2*(cos(theta2)*ihat + sin(theta2)*jhat) ...
        + (l_tip + lc + l_r/2) * (cos(theta3)*ihat+ sin(theta3)*jhat);
    end
end