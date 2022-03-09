function prior = priorUniform(min_vals, max_vals)

% Create a sampling function
prior.sample = @() min_vals + (max_vals - min_vals) .* rand(size(min_vals));

% Create a log prior evaluation function
prior.logpdf = @(theta) log(1 - any(theta < min_vals, 2)) + log(1 - any(theta > max_vals, 2) ) + sum( -log(max_vals - min_vals) );

% Store the bounds also for parameter transform purposes
prior.theta_min = min_vals;
prior.theta_max = max_vals;

end

