clear; close all; 
% Parameters
theta0 = (90- 17/2)*pi/180;
theta_max = pi/4; % Maximum angle
num_frames = 100; % Number of animation frames
theta_vals = linspace(theta0, theta_max, num_frames);

% GIF parameters
gif_filename = 'prbm_animation.gif'; % Output GIF file name
delay_time = 0.1; % Delay between frames (in seconds)

% Example: Two-link manipulator (replace with your robot's kinematics)
lfull = 50; 
ltip = 3.25; 
lrigid = 5; 
ltop = 17; 
lsmall = 1; 

% Find Model Variables from CAD vars
l = lfull - lrigid/2; 
y = 0.8517; 
lc = (1-y)*l/2;

% Initialize figure
figure;
hold on;
axis equal;
grid on;
xlim([-lfull, lfull]);
ylim([0, lfull+ltip]);
xlabel('X-axis');
ylabel('Y-axis');
title('PRBM Animation');

% Initialize plot handles
link1_plot = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % Link 1
link2_plot = plot([0, 0], [0, 0], 'b-', 'LineWidth', 2); % Link 2
link3_plot = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % Link 3

link4_plot = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % Link 1
link5_plot = plot([0, 0], [0, 0], 'b-', 'LineWidth', 2); % Link 2
link6_plot = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % Link 3

link7_plot = plot([0, 0], [0, 0], 'g-', 'LineWidth', 2); % Link 3

joint_plot = plot(0, 0, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k'); % Joints

% Animation loop
for frame = 1:num_frames
    theta = theta_vals(frame); 
    
    theta2 = theta + pi - theta0; 
    theta3 = theta + pi/2 -theta0; 
    ihat = [1 0]; jhat = [0 1];

    rA = (lc+ y*l*sin(theta))*jhat + y*l*cos(theta)*ihat; 
    rtip = rA + ltop/2*(cos(theta2)*ihat + sin(theta2)*jhat) ...
    + (ltip + lc + lrigid/2) * (cos(theta3)*ihat+ sin(theta3)*jhat);
    % Forward kinematics (example for 2-link planar robot)
    % Replace with your own joint/link equations
    r1 =  [lsmall, lc]; 
    r2 = r1+ y*l*sin(theta)*jhat + y*l*cos(theta)*ihat;
    r3 = r2 + (lc + lrigid/2) * (cos(theta3)*ihat+ sin(theta3)*jhat); 
    
    theta = 2*theta0 -theta;

    r4 =  [-lsmall, lc]; 
    r5 = r4+ y*l*sin(theta)*jhat - y*l*cos(theta)*ihat;
    r6 = r5 + (lc + lrigid/2) * (cos(theta3)*ihat+ sin(theta3)*jhat); 

    
    % Update link and joint positions
    set(link1_plot, 'XData', [lsmall, r1(1)], 'YData', [0, r1(2)]);
    set(link2_plot, 'XData', [r1(1), r2(1)], 'YData', [r1(2), r2(2)]);
    set(link3_plot, 'XData', [r2(1), r3(1)], 'YData', [r2(2), r3(2)]);
    set(link4_plot, 'XData', [-lsmall, r4(1)], 'YData', [0, r4(2)]);
    set(link5_plot, 'XData', [r4(1), r5(1)], 'YData', [r4(2), r5(2)]);
    set(link6_plot, 'XData', [r5(1), r6(1)], 'YData', [r5(2), r6(2)]);
    set(link7_plot, 'XData', [r3(1), r6(1)], 'YData', [r3(2), r6(2)]);

    set(joint_plot, 'XData', [r1(1), r2(1), r3(1), r4(1), r5(1), r6(1)],...
        'YData', [r1(2), r2(2), r3(2), r4(2), r5(2), r6(2)]);
   
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
