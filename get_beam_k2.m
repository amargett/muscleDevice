function get_beam_k2(p_in, E, t0, r_well, w, t_min)

    close all; 
    % h : parameters for parametric curve of the inside beam
    R0 = r_well -t0/2; N = 100;  
    rpi = get_r(pi, t_min);
    F = 3E-7; % 300 uN of force   
    
    %% Define theta range
    theta_vals = linspace(0, pi, N);
    
    %% Plot 1: M(theta) for theta_b = pi/2 and theta_b = pi
    figure;
    hold on;
    plot(theta_vals, arrayfun(@(theta) get_M(theta, pi/2), theta_vals), 'b-', 'DisplayName', '\theta_b = \pi/2');
    plot(theta_vals, arrayfun(@(theta) get_M(theta, pi), theta_vals), 'r--', 'DisplayName', '\theta_b = \pi');
    xlabel('\theta [rad]');
    ylabel('M(\theta) [Nm]');
    legend('Location', 'best');
    title('Bending Moment M(\theta)');
    grid on; 
    
    %% Plot 2: Thickness t(theta)
    figure;
    plot(theta_vals, arrayfun(@(theta) get_t(theta), theta_vals), 'k-', 'LineWidth', 1.5);
    xlabel('\theta [rad]');
    ylabel('Thickness t(\theta) [m]');
    title('Thickness vs. \theta');
    grid on;
    

    %% Plot 3: Angular deflection psi(theta_b) from pi/2 to pi
    theta_b_vals = linspace(pi/2, pi, N);
    psi_vals = arrayfun(@(theta_b) get_psi(theta_b), theta_b_vals);

    figure;
    plot(theta_b_vals, psi_vals, 'b-', 'LineWidth', 1.5);
    xlabel('\theta_b [rad]');
    ylabel('\psi(\theta_b) [rad]');
    title('Angular Deflection \psi(\theta_b) from \theta_b = \pi/2 to \pi');
    grid on;

    %% Plot 4: Radial deflection delta(theta_b) from pi/2 to pi
    delta_vals = arrayfun(@(theta_b) get_delta(theta_b), theta_b_vals);

    figure;
    plot(theta_b_vals, delta_vals, 'r-', 'LineWidth', 1.5);
    xlabel('\theta_b [rad]');
    ylabel('\delta(\theta_b) [m]');
    title('Radial Deflection \delta(\theta_b) from \theta_b = \pi/2 to \pi');
    grid on;

    %% Plot 5: Stiffness K from pi/2 to pi
    K_vals = F./delta_vals; 

    figure;
    plot(theta_b_vals, K_vals, 'r-', 'LineWidth', 1.5);
    xlabel('\theta_b [rad]');
    ylabel('K [N/m]');
    title('Stiffness K from \theta_b = \pi/2 to \pi');
    grid on;

 
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
        r_in = 1E-3* get_rp(theta, p_in); 
        
    end

    function r = get_r(theta,t)
        %finds radius at each theta in 
        if theta <= pi/2
            r = R0; 
        else
            r_in = get_r_in(theta); 
            r_out = r_in + t;
            r = (r_in+r_out)./2; 
        end
    end
    
    function [Fbx, Fby]  = get_Fb(F, theta_b)
        Fb = F*(R0 + rpi)/(sin(theta_b)*(rpi+get_r_in(theta_b)*cos(theta_b))...
            -cos(theta_b)*(get_r_in(theta_b)*sin(theta_b))); 
        Fby = Fb*cos(theta_b); Fbx = Fb*sin(theta_b); 
    end


    function Mb = get_Mb(F, theta_b, theta, t)
        [Fbx, Fby] = get_Fb(F, theta_b); 
        r = get_r(theta, t); r_in = get_r_in(theta_b); 

        Mb = Fbx*(r*cos(theta) - r_in*cos(theta_b))...
            + Fby*(r*sin(theta)-r_in*sin(theta_b)); 
    end

    function t = get_t(theta)
        if theta<= pi/2
            t = t0; 
        else
            t = t0 + (t_min - t0) * (theta - pi/2) / pi/2;
        end
    end

    function M = get_M(theta, theta_b)
        t = get_t(theta); 
        r = get_r(theta, t); 
    
        if theta <= theta_b
            M = F*(R0 - r * cos(theta)); 
        else
            M = F*(R0 - r * cos(theta)) + get_Mb(F, theta_b, theta, t); 
        end
    end


    function I = get_I(theta)
        t = get_t(theta);
        I = w*t.^3/12; 
    end

    % Function to compute psi_vals using trapezoidal rule for numerical integration
    function psi_vals = get_psi_vals(theta_b)
        psi_vals = zeros(size(theta_vals)); 
        for i = 2:N
            theta = theta_vals(i); 
            M = get_M(theta, theta_b); 
            I = get_I(theta); 
            
            % Apply trapezoidal rule step (assuming equal spacing in theta_vals)
            h = theta_vals(i) - theta_vals(i-1);  % Step size
            if i == 2
                psi_vals(i) = (M / (E * I)) * h / 2;  % First value has only half weight
            else
                psi_vals(i) = (M / (E * I)) * h;  % Subsequent values have full weight
            end
        end
    end
    
    % Function to compute psi using the trapezoidal rule
    function psi = get_psi(theta_b)
        psi_vals = get_psi_vals(theta_b);
        h = theta_vals(2) - theta_vals(1);  % Step size
        % Apply trapezoidal rule for integration over psi_vals
        psi = (h / 2) * (psi_vals(1) + psi_vals(end) + 2 * sum(psi_vals(2:end-1)));
    end
    
    % Function to compute delta_vals using trapezoidal rule
    function delta_vals = get_delta_vals(theta_b)
        delta_vals = zeros(size(theta_vals)); 
        psi_vals = get_psi_vals(theta_b); 
        for i = 2:N
            h = theta_vals(i) - theta_vals(i-1);  % Step size for delta
            if i == 2
                delta_vals(i) = psi_vals(i) * h / 2;  % First value has half weight
            else
                delta_vals(i) = psi_vals(i) * h;  % Subsequent values have full weight
            end
        end
    end
    
    % Function to compute delta using the trapezoidal rule
    function delta = get_delta(theta_b)
        delta_vals = get_delta_vals(theta_b);
        h = theta_vals(2) - theta_vals(1);  % Step size
        % Apply trapezoidal rule for integration over delta_vals
        delta = (h / 2) * (delta_vals(1) + delta_vals(end) + 2 * sum(delta_vals(2:end-1)));
    end



end