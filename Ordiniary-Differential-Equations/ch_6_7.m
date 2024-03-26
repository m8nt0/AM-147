% Generate synthetic data
x = linspace(0, 2*pi, 100)';
y_true = sin(x) + cos(2*x);  % True underlying function
y_observed = y_true + 0.1*randn(size(x));  % Observed data with noise

% Define basis functions
k = 5;  % Number of basis functions
omega = linspace(0.5, 2.5, k);  % Frequencies of the sinusoidal functions

% Construct design matrix
X_design = zeros(length(x), k);
for j = 1:k
    X_design(:, j) = sin(omega(j) * x);
end

% Solve for coefficients using least squares
coefficients = X_design \ y_observed;

% Generate predicted values using the model
y_predicted = X_design * coefficients;

% Plot the observed data and the fitted model
figure;
plot(x, y_observed, 'bo', 'DisplayName', 'Observed Data');
hold on;
plot(x, y_predicted, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Model');
xlabel('x');
ylabel('y');
title('Regression by Weighted Sum of Sinusoidal Functions');
legend('show');
grid on;
hold off;
