% Program 5.1 Romberg integration
% Computes approximation to definite integral
% Inputs: Matlab function specifying integrand f,
% a, b integration interval, n=number of rows
% Output: Romberg tableau r
function r = romberg(f, a, b, n)
    h = (b - a) ./ (2.^(0:n-1));
    r = zeros(n, n);
    
    % Calculate the first column entries using composite trapezoidal rule
    r(1, 1) = (b - a) * (f(a) + f(b)) / 2;
    
    % Fill in the rest of the tableau using Romberg Integration
    for j = 2:n
        subtotal = 0;
        for i = 1:2^(j-2)
            subtotal = subtotal + f(a + (2*i - 1) * h(j));
        end
        r(j, 1) = r(j-1, 1) / 2 + h(j) * subtotal;
        
        % Perform extrapolation for higher-order approximations
        for k = 2:j
            r(j, k) = (4^(k-1) * r(j, k-1) - r(j-1, k-1)) / (4^(k-1) - 1);
        end
    end
end

% Example 5.11: Apply Romberg Integration to approximate ∫[1 to 2] ln(x) dx
res = romberg(@log, 1, 2, 4);
disp(res);
