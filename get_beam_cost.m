function cost = get_beam_cost(h, params, p_in)
    % h, p_in : parameters for parametric curve of the outside, inside beams
    
    %% Data Tracking Variables
    global data2; 
    %% Get Input Variables
    E = params(1); t0 = params(2); r_well = params(3); w = params(4); 
    R0 = r_well -t0/2;
    theta_f = 3*pi/4; N = 100; 
    %% Section 1
    delta1 = 12*R0^3/(E*w*t0^3)*(3*pi/4-2) ;
    delta2_high = get_delta2(pi/2) ;
    delta2_low = get_delta2(theta_f);
    
    K_high = 1/(delta1+delta2_high);
    K_low = 1/(delta1+delta2_low);

    data2 = [K_high, K_low]; 

    cost = K_low/K_high;

    if isnan(cost)||cost<0
        cost = 10000000; 
    end
    %% Helper Functions
    function radius = get_r(theta, p)
        th = theta-pi/2;

        rootFunc = @(r) (r * cos(th)) - polyval(p, r * sin(th));
        
        % Provide an initial guess for r (e.g., r = 1)
        r_guess = R0*10^3;
        
        % Solve for r using fzero (only finds one root)
        try
            r_sol = fzero(rootFunc, r_guess);
            if r_sol > 0  % Only keep positive radii
                radius = r_sol;
            else 
                radius = NaN;
            end
       
        catch
            radius = NaN;  % No valid solution
        end
    end
    function [r,t] = get_rt(theta)
        r_in = get_r(theta, p_in);
        r_out = get_r(theta, h) ;       
        r = (r_in+r_out)/2;
        t = r_out - r_in;
        if  t<0.05
            t = NaN;% avoids too small thicknesses
        end

    end
    function delta2 = get_delta2(theta_b)
        % gets second section stiffness as a function of boundary condition
        theta_vals1 = linspace(pi/2, theta_b, N);
        theta_vals2 = linspace(theta_b, theta_f, N); 
        r_vals = zeros(size(theta_vals1)); 
        t_vals = zeros(size(theta_vals1)); 
        
        du1 = zeros(size(theta_vals1)); 
        for i = 1:length(theta_vals1)
            theta = theta_vals1(i); 
            [r,t]= get_rt(theta);
            r_vals(i)=r; t_vals(i) = t; 
            du1(i) = (R0-r*cos(theta))^2*r/t^3; 
        end

        du2 = zeros(size(theta_vals2)); 
        for i = 1:length(theta_vals2)
            theta = theta_vals2(i); 
            [r,t]= get_rt(theta);
            du2(i) = (R0-r*cos(theta)-r*sin(theta-theta_b)/ ...
                cos(theta_b -pi/2))^2*r/t^3; 
        end        
        d1 = 12/(E*w)*(trapz(theta_vals1, du1));
        d2 = 12/(E*w)*trapz(theta_vals2, du2);
        delta2 = 12/(E*w)*(trapz(theta_vals1, du1) + trapz(theta_vals2, du2));

        % figure(); plot(theta_vals1, r_vals); title('r'); 
        % figure(); plot(theta_vals1, t_vals); title('t'); 
    end        
end