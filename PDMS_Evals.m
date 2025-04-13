% MATLAB Script to Calculate Young's Modulus from CSV Data
close all
filenames = {'CSVs/1-10_1.csv', 'CSVs/1-10_2.csv'...
    , 'CSVs/1-20_1.csv', 'CSVs/1-20_2.csv', 'CSVs/1-30_1.csv', 'CSVs/1-30_2.csv', 'CSVs/1-30_3.csv'}; 

titles = {'1-10 1', '1-10 2', '1-20 1', '1-20 2', '1-33 1', '1-33 2', '1-33 3'};
% ts = [3.12, 2.98, 1.51, 2.97, 3.05, 2.99, 3.04, 3.02]; 
ts = [3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5]; 
ws = [4.95, 4.95, 4.86, 4.96, 4.98, 5.06, 4.99];
gauge_length = 14.8; % Initial length (mm)

for i = 1:length(filenames)
    get_E(filenames{i}, ts(i)*ws(i), gauge_length, titles{i})
end

%%
function get_E(filename, area, gauge_length, title_text)
    % Load CSV file (update filename if needed)
    opts = detectImportOptions(filename, 'NumHeaderLines', 8); % Skip initial text rows
    data = readtable(filename, opts);
    
    displacement = data{:, 3}; % Displacement (mm)
    force = data{:, 4}; % Force (N)
    
    % Calculate strain 
    strain = (displacement - displacement(1)) / gauge_length;
    
    % Calculate stress
    stress = force / area; % N/mm² = MPa
    
    % Identify the linear region (usually ~0-10% strain)
    linear_region = strain < 0.1; % Modify if needed
    
    % Perform linear fit in the elastic region
    p = polyfit(strain(linear_region), stress(linear_region), 1);
    youngs_modulus = p(1); % Young’s modulus (MPa)
    
    % Display result
    fprintf('Young''s Modulus: %.3f kPa\n', youngs_modulus*10^3);

    disp('filename: ' + string(filename))
    
    % Plot Stress-Strain Curve
    figure;
    plot(strain, stress, 'bo-', 'MarkerSize', 5, 'LineWidth', 1.5);
    hold on;
    plot(strain(linear_region), polyval(p, strain(linear_region)), 'r-', 'LineWidth', 2);
    xlabel('Strain');
    ylabel('Stress (MPa)');
    new_title = title_text + ", E: " + string(youngs_modulus*10^3) + " kPa";
    title(new_title);
    legend('Experimental Data', 'Linear Fit');
    grid on;
end