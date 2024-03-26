function Program_0_1_Nested_Multiplication()
    nested_multiplication();
    decimal_to_binary();
    binary_to_decimal();
    double_precision_representation();
    floating_point_operations();
    loss_of_significance();
    polynomial_evaluation();
    intermediate_value_theorem();
    mean_value_theorem();
    taylor_polynomial();
end

function nested_multiplication()
    % Nested multiplication
    d = 4;
    c = [2, 3, -3, 5, -1];
    x = 0.5;
    b = zeros(d, 1);
    result = nest(d, c, x, b);
    disp(['Nested multiplication result: ' num2str(result)]);
end

function decimal_to_binary()
    % Decimal to binary conversion
    decimalNumber = 25;
    binaryRepresentation = dec2bin(decimalNumber);
    disp(['Binary representation of ' num2str(decimalNumber) ': ' binaryRepresentation]);
end

function binary_to_decimal()
    % Binary to decimal conversion
    binaryNumber = '1101';
    decimalValue = bin2dec(binaryNumber);
    disp(['Decimal value of ' binaryNumber ': ' num2str(decimalValue)]);
end

function double_precision_representation()
    % Double precision representation
    number = 9.4;
    format hex;
    disp(['Hex representation of ' num2str(number) ': ' num2hex(number)]);
    format;
end

function floating_point_operations()
    % Floating-point operations
    a = 1;
    b = 3 * 2^(-53);
    result = (a + b) - 1;
    disp(['Result of (1 + 3 * 2^(-53)) - 1: ' num2str(result)]);
end

function loss_of_significance()
    % Loss of significance
    a = 1.2345678912345678;
    b = 1.2345678912345677;
    result = a - b;
    disp(['Result of (a - b): ' num2str(result)]);
end

function polynomial_evaluation()
    % Polynomial evaluation
    coefficients = [2, 3, -3, 5, -1];
    x_value = 1/2;
    result = polyval(coefficients, x_value);
    disp(['Polynomial evaluation result: ' num2str(result)]);
end

function intermediate_value_theorem()
    % Intermediate Value Theorem
    f = @(x) x.^2 - 3;
    a = 1;
    b = 3;
    y = 0;
    c = fzero(@(x) f(x) - y, [a, b]);
    disp(['Value of c: ' num2str(c)]);
end

function mean_value_theorem()
    % Mean Value Theorem
    f = @(x) x.^2 - 3;
    a = 1;
    b = 3;
    c = fminbnd(@(x) -f(x), a, b);
    slope_c = (f(b) - f(a)) / (b - a);
    disp(['Value of c: ' num2str(c)]);
    disp(['Slope at c: ' num2str(slope_c)]);
end

function taylor_polynomial()
    % Taylor Polynomial
    f = @(x) sin(x);
    k = 4;
    syms x;
    P_k = taylor(f(x), x, 'Order', k);
    x_value = 0.0001;
    remainder_sym = taylor(f(x), x, 'Order', k + 1);
    remainder_func = matlabFunction(remainder_sym);
    remainder_at_x = remainder_func(x_value);
    disp(['Degree ' num2str(k) ' Taylor Polynomial: ' char(P_k)]);
    disp(['Remainder at x = 0.0001: ' num2str(remainder_at_x)]);
end

function y = nest(d, c, x, b)
    if nargin < 4
        b = zeros(d, 1);
    end
    y = c(d + 1);
    for i = d:-1:1
        y = y * (x - b(i)) + c(i);
    end
end
