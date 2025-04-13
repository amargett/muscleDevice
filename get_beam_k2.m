function get_beam_k2(p_in, E, t0, r_well, w, t_min)

    close all; 
    % h : parameters for parametric curve of the inside beam
    R0 = r_well -t0/2; N = 100; step = pi/N; 
    rpi = get_r(pi);
    F = 3E-7; % 300 uN of force   
    
    %% Define theta range
    theta_vals = linspace(0, pi, N);
    theta_b_vals = linspace(pi/2, 99*pi/100, N);
    
    poly_order = 5; 
    p_out = fit_p(); 

    poly_eq = sprintf('%.6f*x^%d', p_out(1), poly_order); % Start with the highest-order term
    for k = 2:length(p_out)
        poly_eq = sprintf('%s + %.8f*x^%d', poly_eq, p_out(k), poly_order - (k - 1));
    end

    disp(poly_eq); 
    
    %% Plot 1: M(theta) for theta_b = pi/2 and theta_b = pi
    figure;
    hold on;
    plot(theta_vals, arrayfun(@(theta) get_M(theta, pi/2), theta_vals), 'b-', 'DisplayName', '\theta_b = \pi/2');
    plot(theta_vals, arrayfun(@(theta) get_M(theta, 9*pi/16), theta_vals), '-', 'DisplayName', '\theta_b = 9\pi/16');
    plot(theta_vals, arrayfun(@(theta) get_M(theta, 5*pi/8), theta_vals), '-', 'DisplayName', '\theta_b = 5\pi/8');
    plot(theta_vals, arrayfun(@(theta) get_M(theta, 3*pi/4), theta_vals), 'g-', 'DisplayName', '\theta_b = 3\pi/4');
    plot(theta_vals, arrayfun(@(theta) get_M(theta, 7*pi/8), theta_vals), '-', 'DisplayName', '\theta_b = 7\pi/8');

    xlabel('\theta [rad]');
    ylabel('M(\theta) [Nm]');
    legend('Location', 'best');
    title('Bending Moment M(\theta)');
    grid on; 

    figure; 
    hold on; 
    plot(theta_b_vals, arrayfun(@(theta_b) get_Fbx(theta_b), theta_b_vals), 'b-', 'DisplayName', 'Fbx');
    plot(theta_b_vals, arrayfun(@(theta_b) get_Fby(theta_b), theta_b_vals), 'r-', 'DisplayName', 'Fby');
    xlabel('\theta_b [rad]');
    ylabel('Fb [N]');
    legend('Location', 'best');
    title('Boundary Condition Forces vs. \theta_b');
    grid on;

    figure; 
    hold on; 
    plot(theta_vals, arrayfun(@(theta) get_r_in(theta), theta_vals), 'b-', 'DisplayName', 'r in');
    plot(theta_vals, arrayfun(@(theta) get_r(theta), theta_vals), '-', 'DisplayName', 'r');
    legend('Location', 'best');
    title('radii');
    grid on
    
    % %% Plot 2: Thickness t(theta)
    % figure;
    % plot(theta_vals, arrayfun(@(theta) get_t(theta), theta_vals), 'k-', 'LineWidth', 1.5);
    % xlabel('\theta [rad]');
    % ylabel('Thickness t(\theta) [m]');
    % title('Thickness vs. \theta');
    % grid on;
    % 
    delta_vals = arrayfun(@(theta_b) get_du_delta(theta_b), theta_b_vals);
    %% Plot 5: Stiffness K from pi/2 to pi
    K_vals = F./delta_vals; 

    figure;
    plot(theta_b_vals, K_vals, 'r-', 'LineWidth', 1.5);
    xlabel('\theta_b [rad]');
    ylabel('K [N/m]');
    title('Stiffness K from \theta_b = \pi/2 to \pi');
    grid on;
    % 
    
    %% Plot 6: Stresses
    figure;
    hold on 
    plot(theta_vals,arrayfun(@(theta) get_stress(theta, pi/2), theta_vals), 'b-', 'DisplayName', '\theta_b = \pi/2');
    plot(theta_vals,arrayfun(@(theta) get_stress(theta, pi), theta_vals), 'r-', 'DisplayName', '\theta_b = \pi');
    xlabel('\theta [rad]');
    ylabel('stress \sigma) [N/m^2]');
    legend('Location', 'best');
    title('Stress');
    grid on;
 
    %% Helper Functions
    function stress = get_stress(theta, theta_b)
        M = get_M(theta, theta_b); 
        t = get_t(theta);
        c = t/2; 
        I = w*t^3/12; 

        stress= M*c/I; 
    end

    function p_out = fit_p()
        x_out = zeros(1, N/2); y_out = zeros(1, N/2);
        for i = 1:N/2
            th = theta_vals(i+N/2); 
            t = get_t(th); r = get_r_in(th) + t; 
            x_out(i) = 10^3*r*sin(th-pi/2); y_out(i) = 10^3*r*cos(th-pi/2); 
        end

        p_out = polyfit(x_out, y_out, poly_order);

        % Evaluate the polynomial for a smooth curve
        x_fit = linspace(min(x_out), max(x_out), 100);
        y_fit = polyval(p_out, x_fit);
        y_fit2 = polyval(p_in, x_fit); 
                
        % Plot original data
        figure;
        plot(x_out, y_out, 'o', 'DisplayName', 'Outer Boundary Original Data');
        hold on;
        
        % Plot polynomial fit
        plot(x_fit, y_fit, '-', 'DisplayName', ['Outer Boundary Poly Fit (Order ', num2str(poly_order), ')']);

        % Plot polynomial fit
        plot(x_fit, y_fit2, '-', 'DisplayName', ['PRBM Poly Fit (Order ', num2str(poly_order), ')']);
        
        legend('Location','best'); 
        xlabel('X position [mm]'); ylabel('Y position [mm]')
        hold off
    end

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

    function Fbx = get_Fbx(theta_b)
        [Fbx, ~, ~] = get_Fb(theta_b); 
    end
    function Fby = get_Fby(theta_b)
        [~, Fby, ~] = get_Fb(theta_b);
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
       if theta_b == 99*pi/100
           F/delta
       end
    end
end