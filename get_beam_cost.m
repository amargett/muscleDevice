function cost = get_beam_cost(h, params, p_in, last)
    % h, p_in : parameters for parametric curve of the outside, inside beams
    
    %% Data Tracking Variables
    global data2; 
    %% Get Input Variables
    E = params(1); t0 = params(2); r_well = params(3); w = params(4); t_min = params(5); 
    R0 = r_well -t0/2; N = 100; [rpi, tpi] = get_rt(pi); 

    if last 
        disp(rpi); disp(tpi); 
    end
    %% Find Stiffness
    K_high = get_K(pi/2); K_low = get_K(pi); data2 = [K_high, K_low]; 

    cost = K_low/K_high; % cost function

    if isnan(cost)||cost<0
        cost = 10000000; 
    end
    %% Helper Functions
    function radius = get_r(theta, p)
        % finds radius from a polynomial function
        th = theta-pi/2;

        rootFunc = @(r) (r * cos(th)) - polyval(p, r * sin(th));

        try
            radius = fzero(rootFunc, R0*10^3); % Solve for r using fzero
        catch
            radius = NaN;  % No valid solution
            if last 
                disp('no valid sol'); 
            end
        end
        if radius < 0 
            radius = NaN; % only keeps positive radii
            if last
                disp('negative radius');
            end
        end
    end

    function [r,t] = get_rt(theta)
        %finds radius and thickness at each theta in m
        if theta < pi/2
            r = R0; 
            t = t0; 
        else
            r_in = get_r(theta, p_in); r_out = get_r(theta, h); 
            % find r, t, convert from mm to m
            r = 1E-3*(r_in+r_out)/2; t = 1E-3*(r_out - r_in);
            if  t<t_min
                if last
                    disp('thickness too small'); 
                    disp(t); disp(r_in); disp(r_out); 
                end
                t = NaN;% avoids too small thicknesses
            end
        end
    end

    function K = get_K(theta_b)
        % gets beam stiffness as a function of boundary condition'
        theta_vals = linspace(0, pi, N); du = zeros(size(theta_vals)); 
        t_vals = zeros(size(theta_vals)); r_vals = zeros(size(theta_vals));

        for i = 1:length(theta_vals) 
            theta = theta_vals(i); [r,t] = get_rt(theta);
            t_vals(i) = t; r_vals(i)= r; 
            if theta < theta_b
                M = R0 -r*cos(theta); 
            else
                M = R0 -r*cos(theta)+ r*cos(theta)*(rpi+R0)/rpi; 
            end
            du(i) = M^2*r/t^3;     
        end
        
        K = w*E/(12*trapz(theta_vals, du)); 
        if last
            plot(theta_vals, t_vals); title('t'); 
            figure(); plot(theta_vals, r_vals); title('r'); 
        end
    end   
end