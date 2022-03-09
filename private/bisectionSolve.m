function x_mid = bisectionSolve( f, window )
% Attempts to find a zero to input function f, within the window provided.

max_iters = 1000;
tol = 1e-6;

% Initialise the method
x_left = window(1);
x_right = window(2);
f_left = f(x_left);
f_right = f(x_right);
% Check for a suitable window
if sign(f_left) * sign(f_right) > 0
    error('Bisection method window does not contain a zero!');
end

% Loop until maximum iterations or 
looping = true; iters = 0;
while looping

    % Increment iteration count
    iters = iters + 1;
    
    % Find and evaluate midpoint
    x_mid = (x_left + x_right) / 2;
    f_mid = f(x_mid);
    
    % Adjust window to use the midpoint so it still contains a zero
    if f_left * f_mid < 0
        x_right = x_mid;
        f_right = f_mid;
    else
        x_left = x_mid;
        f_left = f_mid;
    end
    
    % Check loop criteria
    if iters == max_iters || abs(f_mid) < tol
        looping = false;
    end
        
end

if iters == max_iters
    warning('Bisection method failed to converge to specified tolerance!');
end