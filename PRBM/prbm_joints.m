clear; close all; 
% Parameters
theta_max = pi/4; % Maximum angle
num_frames = 100; % Number of animation frames
theta_vals = linspace(0, pi/2, num_frames);

% GIF parameters
gif_filename = 'prbm_animation.gif'; % Output GIF file name
delay_time = 0.05; % Delay between frames (in seconds)

%% Parameters to set
r_well = 7.5; % mm
l = 3.75; %mm
theta = 3*pi/8;
t0 = 1; % mm, thickness of constant part of beam
r_inner = 0.5; % mm, radius of inner circle
l_r = r_well-r_inner-t0-l- 0.5; % rigid length
l_tip = r_well-t0-(l+l_r+r_inner)*sin(theta); % mm, length of variable stiffness tip

%% PRBM/geometric parameters
gamma = 0.8517; 
l_eff = l + l_r/2;  % distance to force application (PRBM)
lc=(1-gamma)*l_eff/2;
l_top = (r_inner+l+l_r)*cos(theta);  


% Initialize figure
figure;
hold on;
axis equal;
grid on;
xlim([-r_well, r_well]);
ylim([-r_well, r_well]);
xlabel('X-axis');
ylabel('Y-axis');
title('PRBM Animation');

% Initialize plot handles
link1_plot = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % Link 1
link2_plot = plot([0, 0], [0, 0], 'b-', 'LineWidth', 2); % Link 2
link3_plot = plot([0, 0], [0, 0], 'g-', 'LineWidth', 2); % Link 3
link4_plot = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % Link 4
link5_plot = plot([0, 0], [0, 0], 'b-', 'LineWidth', 2); % Link 5

joint_plot = plot(0, 0, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k'); % Joints

% Animation loop
for frame = 1:num_frames
    theta_p = theta_vals(frame);

    ihat = [1 0]; jhat = [0 1];

    % Forward kinematics (example for 2-link planar robot)
    % Replace with your own joint/link equations
    r1 = (r_inner +lc)*(cos(theta)*ihat + sin(theta)*jhat);
    r2 = r1 + gamma*l*(cos(theta-theta_p)*ihat +sin(theta-theta_p)*jhat);
    r3 = r2 + (lc+l_r)*(cos(theta)*ihat + sin(theta)*jhat); 
    r4 = r3 + l_top*(cos(theta_p)*-ihat + sin(theta_p)*jhat); 
    r5 = r4 +(l_tip)*(sin(theta_p)*ihat + cos(theta_p)*jhat);

    % Update link and joint positions
    set(link1_plot, 'XData', [0, r1(1)], 'YData', [0, r1(2)]);
    set(link2_plot, 'XData', [r1(1), r2(1)], 'YData', [r1(2), r2(2)]);
    set(link3_plot, 'XData', [r2(1), r3(1)], 'YData', [r2(2), r3(2)]);
    set(link4_plot, 'XData', [r3(1), r4(1)], 'YData', [r3(2), r4(2)]);
    set(link5_plot, 'XData', [r4(1), r5(1)], 'YData', [r4(2), r5(2)]);

    set(joint_plot, 'XData', [r1(1), r2(1), r3(1), r4(1), r5(1)],...
                    'YData', [r1(2), r2(2), r3(2), r4(2), r5(2)]);
   
    % Capture the frame
    frame_image = getframe(gcf);
    im = frame2im(frame_image); % Convert to image
    
    % Write to GIF
    [imind, cm] = rgb2ind(im, 256); % Convert to indexed image
    if frame == 1
        % Create the GIF file on the first frame
        imwrite(imind, cm, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay_time);
    else
        % Append to the GIF file
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
    
    % Pause for real-time animation (optional)
    pause(delay_time);
end

disp(['Animation saved as ', gif_filename]);
