function K = get_beam_k(p_in, E, t0, r_well, w, t_min, theta_b)
    % h : parameters for parametric curve of the inside beam
    %% Get Input Variables
    R0 = r_well -t0/2; N = 100;  
    %% Find Stiffness
    rpi = get_r(pi, t_min);
    K = get_K(theta_b); 

    %% Helper Functions
    function radius = get_rp(theta, p)
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

    function r = get_r(theta,t)
        %finds radius at each theta in 
        r_in =1E-3* get_rp(theta, p_in); r_out = r_in + t; 
        r = (r_in+r_out)/2; 
    end

    function K = get_K(theta_b)
        % gets beam stiffness as a function of boundary condition'
        theta_vals = linspace(0, pi, N); du = zeros(size(theta_vals)); 
        t_vals = linspace(t0, t_min, N/2); r_vals = zeros(size(theta_vals));

        for i = 1:length(theta_vals) 
            theta = theta_vals(i);
            if i < N/2 +1
                r = R0; t = t0; 
            else
            t = t_vals(i-N/2); 
            r = get_r(theta, t);
            end
            r_vals(i)= r; 
            if theta < theta_b
                M = R0 -r*cos(theta); 
            else
                M = R0 -r*cos(theta)+ r*cos(theta)*(rpi+R0)/rpi; 
            end
            du(i) = M^2*r/t^3;     
        end
        
        K = w*E/(12*trapz(theta_vals, du)); 
    end   
end