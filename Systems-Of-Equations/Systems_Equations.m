%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = backSubstitution(L, U, P, b)
    [m, n] = size(U);
    y = zeros(m, 1);
    x = zeros(n, 1);

    % Forward substitution (Ly = Pb)
    for i = 1:m
        y(i) = b(P(i));
        for j = 1:i-1
            y(i) = y(i) - L(i, j) * y(j);
        end
    end

    % Back substitution (Ux = y)
    for i = m:-1:1
        x(i) = y(i);
        for j = i+1:n
            x(i) = x(i) - U(i, j) * x(j);
        end
        x(i) = x(i) / U(i, i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given function and its fixed points
g = @(x) sqrt(6./(1+x));
fixed_points = fzero(@(x) g(x) - x, [0, 2]);

disp(['Fixed Points: ' num2str(fixed_points)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualize the geometric interpretation
x_vals = linspace(0, 2, 100);
y_vals = g(x_vals);

figure;
plot(x_vals, y_vals, 'b', 'LineWidth', 2);
hold on;
plot(x_vals, x_vals, 'r--', 'LineWidth', 2);
title('Geometry of Fixed-Point Iteration');
xlabel('x_{n}');
ylabel('x_{n+1}');
legend('g(x)', 'y=x', 'Location', 'Best');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial guess
x0 = 1;

% Apply Fixed-Point Iteration
[x, convergence_rate] = fixed_point_iteration(g, x0, tolerance, max_iterations);

disp(['Converged Root: ' num2str(x)]);
disp(['Convergence Rate: ' num2str(convergence_rate)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 1.2 Fixed-Point Iteration
% Computes approximate solution of g(x) = x
% Input: function handle g, starting guess x0, number of iteration steps k
% Output: Approximate solution xc

function xc = fpi(g, x0, k)
    x(1) = x0;
    for i = 1:k
        x(i+1) = g(x(i));
    end
    xc = x(k+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the function g(x) = cos(x)
g = @(x) cos(x);

% Call the Fixed-Point Iteration function
xc = fpi(g, 0, 10);

% Display the result
disp(['Approximate solution: ', num2str(xc)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial guess
x0 = 1;

% Apply Fixed-Point Iteration
[x, convergence_rate] = fixed_point_iteration(g, x0, tolerance, max_iterations);

disp(['Converged Root: ' num2str(x)]);
disp(['Convergence Rate: ' num2str(convergence_rate)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to find the root of
f = @(x) x.^3 - 6*x.^2 + 11*x - 6;

% Bracketing interval [a, b]
a = 0;
b = 3;

% Plot the function to visualize the bracketing interval
x_vals = linspace(a, b, 100);
figure;
plot(x_vals, f(x_vals));
hold on;
plot([a, b], [0, 0], 'r', 'LineWidth', 2);
title('Bracketing a Root');
xlabel('x');
ylabel('f(x)');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tolerance for accuracy
tolerance = 1e-6;

% Maximum number of iterations
max_iterations = 100;

% Bisection method
[root, iterations] = bisection(f, a, b, tolerance, max_iterations);

disp(['Root: ' num2str(root)]);
disp(['Iterations: ' num2str(iterations)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bisection Method Demo

% Define the function
f = @(x) x^3 + x - 1;

% Initial interval [a, b]
a = 0;
b = 1;

% Tolerance
tol = 1e-6;

% Bisection Method
xc = bisect(f, a, b, tol);

% Display the result
fprintf('Approximate root: %.6f\n', xc);

% Bisection Method function
function xc = bisect(f, a, b, tol)
    if sign(f(a)) * sign(f(b)) >= 0
        error('f(a)f(b)<0 not satisfied!');
    end

    fa = f(a);
    fb = f(b);

    while (b - a) / 2 > tol
        c = (a + b) / 2;
        fc = f(c);

        if fc == 0
            break; % c is a solution, done
        end

        if sign(fc) * sign(fa) < 0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        end
    end
end

xc = (a + b) / 2; % new midpoint is the best estimate
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given function and its derivative
f = @(x) x^3 - 6x^2 + 11x - 6;
df = @(x) 3x^2 - 12x + 11;

% Initial guess
x0 = 2;

% Apply Newton's Method
[root, iterations] = newtons_method(f, df, x0, tolerance, max_iterations);

disp(['Root: ' num2str(root)]);
disp(['Iterations: ' num2str(iterations)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given function and its derivative
f = @(x) exp(x) - 2;
df = @(x) exp(x);

% Initial guess
x0 = 0;

% Apply Newton's Method
[root, iterations] = newtons_method(f, df, x0, tolerance, max_iterations);

disp(['Root: ' num2str(root)]);
disp(['Iterations: ' num2str(iterations)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [root, iterates, errors] = newtonsMethod(initialGuess, maxIterations, tolerance)
% Define the function and its derivatives
syms x;
f = x^3 + x - 1;
f_prime = diff(f, x);
f_double_prime = diff(f_prime, x);

scss
Copy code
% Initialize variables
root = initialGuess;
iterates = [root];
errors = [];

% Perform Newton's Method iterations
for i = 1:maxIterations
    % Check for convergence
    if abs(subs(f, x, root)) < tolerance
        break;
    end

    % Update the root using Newton's Method formula
    root = root - subs(f, x, root) / subs(f_prime, x, root);

    % Store the current iterate and error
    iterates = [iterates, root];
    errors = [errors, abs(root - subs(f, x, root))];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialGuess = -0.7;
maxIterations = 100;
tolerance = 1e-8;

[root, iterates, errors] = newtonsMethod(initialGuess, maxIterations, tolerance);

% Display results
disp('Root:');
disp(root);
disp('Iterates:');
disp(iterates);
disp('Errors:');
disp(errors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, b] = naiveGaussianElimination(A, b)
[m, n] = size(A);

css
Copy code
for k = 1:min(m-1, n)
    for i = k+1:m
        factor = A(i, k) / A(k, k);
        A(i, k:n) = A(i, k:n) - factor * A(k, k:n);
        b(i) = b(i) - factor * b(k);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function opCount = operationCounts(A)
[m, n] = size(A);
opCount = 0;

ruby
Copy code
for k = 1:min(m-1, n)
    for i = k+1:m
        opCount = opCount + 2 * n;  % Subtracting and multiplying
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, U, P] = gaussianElimination(A)
[m, n] = size(A);
L = eye(m);
U = A;
P = eye(m);

scss
Copy code
for k = 1:min(m-1, n)
    [~, pivot] = max(abs(U(k:m, k)));
    pivot = pivot + k - 1;

    if pivot ~= k
        % Swap rows in U
        temp = U(k, :);
        U(k, :) = U(pivot, :);
        U(pivot, :) = temp;

        % Swap rows in L
        temp = L(k, 1:k-1);
        L(k, 1:k-1) = L(pivot, 1:k-1);
        L(pivot, 1:k-1) = temp;

        % Swap rows in P
        temp = P(k, :);
        P(k, :) = P(pivot, :);
        P(pivot, :) = temp;
    end

    for i = k+1:m
        factor = U(i, k) / U(k, k);
        L(i, k) = factor;
        U(i, k:n) = U(i, k:n) - factor * U(k, k:n);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = backSubstitution(L, U, P, b)
[m, n] = size(U);
y = zeros(m, 1);
x = zeros(n, 1);

css
Copy code
% Forward substitution (Ly = Pb)
for i = 1:m
    y(i) = b(P(i));
    for j = 1:i-1
        y(i) = y(i) - L(i, j) * y(j);
    end
end

% Back substitution (Ux = y)
for i = m:-1:1
    x(i) = y(i);
    for j = i+1:n
        x(i) = x(i) - U(i, j) * x(j);
    end
    x(i) = x(i) / U(i, i);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example usage
A = [1 2 -1; 2 1 -2; -3 1 1];
b = [3; 3; -6];

[L, U, P] = gaussianElimination(A);
x = backSubstitution(L, U, P, b);
disp(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, U] = matrixGaussianElimination(A)
[m, n] = size(A);
L = eye(m);
