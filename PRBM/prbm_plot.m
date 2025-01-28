close all; clear; 
N = 100; 
theta_vals = linspace((90- 17/2)*pi/180, 0, N); 

rtipdata = zeros(2,N); 
for i = 1:N
    rtipdata(:,i) = prbm(theta_vals(i)); 
end

% Extract x and y components
x_vals = rtipdata(1, :);
y_vals = rtipdata(2, :);

% Fit a polynomial to the data
poly_order = 4; % Change the order of the polynomial if needed
p = polyfit(x_vals, y_vals, poly_order);

% Evaluate the polynomial for a smooth curve
x_fit = linspace(min(x_vals), max(x_vals), 100);
y_fit = polyval(p, x_fit);

% Plot original data
figure;
plot(x_vals, y_vals, 'o', 'DisplayName', 'Original Data');
hold on;

% Plot polynomial fit
plot(x_fit, y_fit, '-', 'DisplayName', ['Polynomial Fit (Order ', num2str(poly_order), ')']);
hold off;

% Add labels and legend
xlabel('X Values');
ylabel('Y Values');
title('Data and Polynomial Fit');
legend('Location', 'Best');
grid on;

poly_eq = sprintf('y = %.6fx^%d', p(1), poly_order); % Start with the highest-order term
for k = 2:length(p)
    poly_eq = sprintf('%s + %.8f*x^%d', poly_eq, p(k), poly_order - (k - 1));
end

disp(poly_eq)