function K_vals = get_beam_kvals(p_in, E, t0, r_well, w, t_min)

    % close all; 
    % h : parameters for parametric curve of the inside beam
    R0 = r_well -t0/2; N = 100;
    rpi = get_r(pi);
    F = 3E-7; % 300 uN of force   
    
    %% Define theta range
    theta_vals = linspace(0, pi, N);
    theta_b_vals = linspace(pi/2, 99*pi/100, N); 

    delta_vals = arrayfun(@(theta_b) get_du_delta(theta_b), theta_b_vals);
    %% Plot 5: Stiffness K from pi/2 to pi
    K_vals = F./delta_vals; 

    %% Helper Functions

    function radius = get_rp(theta, p)
        % finds radius from a polynomial function
        th = theta-pi/2;

        rootFunc = @(r) (r * cos(th)) - polyval(p, r * sin(th));

        try
            radius = fzero(rootFunc, R0*10^3); % Solve for r using fzero
        catch
            disp('No valid solution !')
            radius = 0;  % No valid solution
        end
        if radius < 0 
            disp('radius negative !')
            radius = 0; % only keeps positive radii
        end
    end

    function r_in = get_r_in(theta)
        if theta< pi/2
            r_in = R0 - t0/2; 
        else
            r_in = 1E-3* get_rp(theta, p_in); 
        end
    end

    function r = get_r(theta)
        %finds radius at each theta in 
        t = get_t(theta); 
        if theta <= pi/2
            r = R0; 
        else
            r_in = get_r_in(theta); 
            r_out = r_in + t;
            r = (r_in+r_out)./2; 
        end
    end
    
    function [Fbx, Fby, M_R]  = get_Fb(theta_b)
        r_in = get_r_in(theta_b); 
        r_BR_x = r_in*sin(theta_b); r_BR_y = rpi+ r_in*cos(theta_b); 
        r_FB_y  = R0 - r_in*cos(theta_b); 

        Fb = F*(r_FB_y + r_BR_y)/(sin(theta_b)*r_BR_y - cos(theta_b)*r_BR_x);
        Fby = Fb*cos(theta_b);
        Fbx = Fb*sin(theta_b);

        M_R = Fbx*r_BR_y - F*(R0 + rpi) - Fby*r_BR_x; 
    end

    function M_r = get_MR(theta_b)
        [~, ~, M_r] = get_Fb(theta_b);
    end

    function Mb = get_Mb(theta_b, theta)
        [Fbx, Fby, ~] = get_Fb(theta_b); 
        r = get_r(theta); r_in = get_r(theta_b); 
        r_Bth_y = -r*cos(theta) + r_in*cos(theta_b); 
        r_Bth_x = -r*sin(theta) + r_in*(sin(theta_b)); 

        Mb = -Fbx*r_Bth_y + Fby*r_Bth_x; 
    end

    function t = get_t(theta)
        if theta<= pi/2
            t = t0; 
        else
            t = t0 + (t_min - t0) * (theta - pi/2) / (pi/2);
        end
    end

    function M = get_M(theta, theta_b)
        r = get_r(theta); 
    
        if theta <= theta_b
            M = F*(R0 - r * cos(theta)); 
        else
            M = F*(R0 - r * cos(theta)) + get_Mb(theta_b, theta); 
        end

        if theta == pi
            M = F*(R0 - r * cos(theta)) + get_Mb(theta_b, theta) ...
                + get_MR(theta_b); 
        end
    end

    function I = get_I(theta)
        t = get_t(theta);
        I = w*t^3/12; 
    end
    
    function delta = get_du_delta(theta_b)
        du = zeros(size(theta_vals)); 
        for i = 1:N
            theta = theta_vals(i); 
            M = get_M(theta, theta_b); 
            r = get_r(theta); 
            I = get_I(theta); 
            du(i) = M^2*r/(E*I); 
        end
       U = trapz(theta_vals, du); 
       delta = U/F; 
    end
end