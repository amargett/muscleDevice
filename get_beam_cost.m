function cost = get_beam_cost(h, params, p_in)
    % h, p_in : parameters for parametric curve of the outside, inside beams

    % Initiliazing global variables to track Kb, Kr, p
    
    %% Get Input Variables
    E = params(1); t0 = params(2); r_well = params(3); w = params(4); 
    R0 = r_well +t0/2; 
    theta_f = 7*pi/8; N = 100; 
    
    %% Section 1
    k1 = 12*R0^3/E*w*t0^3*(3*pi/4-2);
    k2_high = get_k2(pi/2);
    k2_low = get_k2(theta_f);

    cost = (k1+k2_low)/(k1+k2_high) ;
    if isnan(cost)
        cost = 10000000; 
    end
    %% Helper Functions
    function r = get_r(theta, p)
        theta_c = theta-pi/2;
        % f = @(r) r*cos(theta_c) - polyval(p, r*sin(theta_c)); % Define function in r
        % options = optimset('Display', 'off'); % Turn off display for fzero
        r = fminbnd(@(r) abs(r*cos(theta_c) - polyval(p, r*sin(theta_c))), 1E-6, 10);
    end
    function [r,t] = get_rt(theta)
        r_in = get_r(theta, p_in);
        r_out = get_r(theta, h);
        r = (r_in+r_out)/2; 
        t = r_out - r_in; 
    end
    function k2 = get_k2(theta_b)
        % gets second section stiffness as a function of boundary condition
        theta_vals1 = linspace(pi/2, theta_b, N);
        theta_vals2 = linspace(theta_b, theta_f, N); 
        du1 = zeros(size(theta_vals1)); 
        for i = 1:length(theta_vals1)
            theta = theta_vals1(i); 
            [r,t]= get_rt(theta);
            du1(i) = (R0-r*cos(theta))^2*r/t^3; 
        end

        du2 = zeros(size(theta_vals2)); 
        for i = 1:length(theta_vals2)
            theta = theta_vals2(i); 
            [r,t]= get_rt(theta);
            du2(i) = (R0-r*cos(theta)-r*sin(theta-theta_b)/ ...
                cos(theta_b -pi/2))^2*r/t^3; 
        end        
        % du1
        % du2
        k2 = 12/(E*w)*(trapz(theta_vals1, du1) + trapz(theta_vals2, du2)) ;
    end        
end