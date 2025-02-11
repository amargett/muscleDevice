function cost = get_prbm_cost(z, params)
    % Initiliazing global variables to track Kb, Kr, p
    global data
    global p
    global x_vals
    global y_vals
    
    %% Get Input Variables
    E = params(1); r_well = params(2); t0 = params(3); r_inner = params(4); 
    poly_order = params(5); SF = params(6); Tmax = params(7); sigma_yield = params(8); 

    theta = z(1); l = z(2); t = z(3); w = z(4); 
    
    %% PRBM/geometric parameters
    gamma = 0.8517; 
    l_r = r_well-r_inner-t0-l; % rigid length
    l_tip = r_well-t0-(l+l_r+r_inner)*sin(theta); % mm, length of variable stiffness tip
    l_eff = l + l_r/2;  % distance to force application (PRBM)
    lc=(1-gamma)*l_eff/2;
    l_top = (r_inner+l+l_r)*cos(theta); 

    %% Fit Polynomial to PRBM Data
    N = 100; theta_vals = linspace(0, pi/2, N); 
    
    rtipdata = zeros(2,N); 
    for i = 1:N
        rtipdata(:,i) = prbm(theta_vals(i)); 
    end
    
    % Extract x and y components
    x_vals = rtipdata(1, :)*10^3; y_vals = rtipdata(2, :)*10^3;

    % Fit a polynomial to the data
    p = polyfit(x_vals, y_vals, poly_order);
    
    %% Radial Stiffness Kr
    I = @(t,w) w*t^3/12; 
    Kr = 1/(sin(theta)^2*l/(2*E*I(t,w))+...
        cos(theta)^2*l^3/(6*E*I(t,w))); 
    
    %% Bending Stiffness Kb
    K_theta = 2.65; alpha = 0.8517; % PRBM coefficients
    Kb = 2*alpha*K_theta*E*I(t,w)/(l_eff); 

    % Kb Upper Bound due to Motor Specs
    Kb_max_motor = 2*Tmax/(SF*pi/2); 
    
    % Kb Upper Bound due to yield stress
    s = get_s();  
    Kb_max_yield = sigma_yield*w*t^2/(SF*6*l_eff*s); 
    
    % penalizes heavily for breaking bounds
    penalty = 100 * (Kb > Kb_max_motor | Kb > Kb_max_yield);
    
    cost = Kb/Kr + penalty; 
    data = [Kb Kr Kb_max_motor Kb_max_yield s];

    %% Helper Functions
    function rC = prbm(theta_p)
        % finds endpoint position as a function of theta
        ihat = [1 0]; jhat = [0 1]; 
        rA = (r_inner +lc)*(cos(theta)*ihat + sin(theta)*jhat); 
        rB = rA + gamma*l_eff*(cos(theta-theta_p)*ihat +sin(theta-theta_p)*jhat); 
        rC = rB + (lc+l_r/2)*(cos(theta)*ihat + sin(theta)*jhat)...
             +(l_tip)*(sin(theta_p)*ihat + cos(theta_p)*jhat)...
            + l_top*(cos(theta_p)*-ihat + sin(theta_p)*jhat);
    end

    function s = get_s()
        s = 0; 
        for j = 1:length(theta_vals)-1
            r = sqrt(x_vals(j)^2 + y_vals(j)^2); 
            s = s + r*(theta_vals(j+1)-theta_vals(j)); 
        end
    end
end