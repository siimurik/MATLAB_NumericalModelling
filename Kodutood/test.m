% Fixed-Point Iteration with relaxation for solving the system of equations
% y + 5xy = x^2
% y + x^2 - x = 0.75

% Define tolerance, maximum number of iterations, and relaxation factor
tolerance = 1e-6;
max_iter = 100;
lambda = 0.5;  % Relaxation factor (0 < lambda <= 1)

% Initial guesses for x and y
x = 0.5;
y = 0.5;

% Fixed-point iteration loop
for iter = 1:max_iter
    % Save previous values of x and y for comparison
    x_prev = x;
    y_prev = y;
    
    % Update y using the second equation in fixed-point form
    y_new = 0.75 - x_prev^2 + x_prev;
    
    % Check if y is a valid number
    if isnan(y_new) || isinf(y_new)
        fprintf('Invalid value encountered for y at iteration %d.\n', iter);
        break;
    end
    
    % Compute the value inside the square root for x
    sqrt_arg = y_new + 5 * x_prev * y_new;
    
    % Ensure we are not taking the square root of a negative number
    if sqrt_arg < 0
        fprintf('Encountered negative value under square root at iteration %d.\n', iter);
        break;
    end
    
    % Update x using the first equation in fixed-point form
    x_new = sqrt(sqrt_arg);
    
    % Check if x is a valid number
    if isnan(x_new) || isinf(x_new)
        fprintf('Invalid value encountered for x at iteration %d.\n', iter);
        break;
    end
    
    % Apply relaxation
    x = (1 - lambda) * x_prev + lambda * x_new;
    y = (1 - lambda) * y_prev + lambda * y_new;
    
    % Check for convergence
    if abs(x - x_prev) < tolerance && abs(y - y_prev) < tolerance
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
    
    % Display iteration results
    fprintf('Iteration %d: x = %.6f, y = %.6f\n', iter, x, y);
end

if iter == max_iter
    fprintf('Max iterations reached without convergence.\n');
end

% Display the final solution
fprintf('Solution: x = %.6f, y = %.6f\n', x, y);
