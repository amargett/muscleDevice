function plot_kvals(p, E, t0, r_well, z)
    %% Plot 5: Stiffness K from pi/2 to pi
    
    t_min_vals = [0.05 0.75 0.1]; % mm
    w_vals = [1 1.5 2 2.5 3.2]; 
    theta_b_vals = linspace(pi/2, pi, 100); 
    figure; 
   
    hold on; % Ensures multiple plots appear on the same figure

    for w = w_vals
        K_vals = get_beam_kvals(p, E, t0, r_well, w*10^-3, 0.1*10^-3); 
        plot(theta_b_vals, K_vals, 'DisplayName', sprintf('w = %.1f mm', w));
    end
    xlabel('\theta_b [rad]');
    ylabel('K [N/m]');
    title('Stiffness K from \theta_b = \pi/2 to \pi');
    grid on;
    legend('Location', 'best');
end
