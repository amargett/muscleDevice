function cost = get_prbm_cost(z, params)
    % Initiliazing global variables to track Kb, Kr, p
    global data
    global p
    data = []; p = []; 

    %% Get Input Variables
    E = params(1); r_well = params(2); t0 = params(3); r_inner = params(4); 
    l_tip = params(5); SF = params(6); Tmax = params(7); sigma_yield = params(8); 

    theta = z(1); l = z(2); t = z(3); w = z(4); 
    
    %% PRBM/geometric parameters
    gamma = 0.8517; 
    l_r = r_well - r_inner - t0 - l_tip -l; % rigid length
    l_eff = l + l_r/2;  % distance to force application (PRBM)
    lc=(1-gamma)*l_eff;
    l_top = (r_inner+l+l_r)*cos(theta); 

    %% Fit Polynomial to PRBM Data
    N = 100; theta_vals = linspace(0, pi/2, N); 
    
    rtipdata = zeros(2,N); 
    for i = 1:N
        rtipdata(:,i) = prbm(theta_vals(i)); 
    end
    
    % Extract x and y components
    x_vals = rtipdata(1, :); y_vals = rtipdata(2, :);

    % Fit a polynomial to the data
    poly_order = 4; % Change the order of the polynomial if needed
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
    this_xy = get_xy(pi/2); this_x = this_xy(1); 
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
        rB = rA + gamma*l_eff*(cos(theta+theta_p)*ihat +sin(theta+theta_p)*jhat); 
        rC = rB + (lc+l_r/2+l_tip)*(cos(theta)*ihat + sin(theta)*jhat)...
            + l_top*(cos(theta+pi/2)*ihat + sin(theta+pi/2)*jhat);
    end

    function xy = get_xy(theta_c)
        % Finds xy array for a given polynomial function at a given theta
        
        eqns = @(v)[ v(1) * tan(theta_c - pi/2) - v(2);  % First equation
                     polyval(p, v(1)) - v(2) ];         % Second equation
        
        initial_guess = [0.5, 0.5];  % Initial guess for [x, y]
        
        % Call fsolve to solve the system
        options = optimset('Display', 'off');  % Suppress display of output
        [xy, ~] = fsolve(eqns, initial_guess, options);
        
    end

    function r = get_r(theta_c)
        xy = get_xy(theta_c); x = xy(1); y = xy(2); 
        r = sqrt(x^2 + y^2); 
    end
    
    function s = get_s()
        r_vals = arrayfun(@(t) get_r(t), theta_vals); % Evaluate r at each theta
        s = 0; 
        for j = 1:length(theta_vals)-1
            s = s + r_vals(j)*(theta_vals(j+1)-theta_vals(j)); 
        end
    end

end