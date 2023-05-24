function dist = distUniform(min_vals, max_vals)

% Create a sampling function
dist.sample = @() min_vals + (max_vals - min_vals) .* rand(size(min_vals));

% Create a log pdf evaluation function
dist.logpdf = @(x) log(1 - any(x < min_vals, 2)) + log(1 - any(x > max_vals, 2) ) + sum( -log(max_vals - min_vals) );

% Store the bounds also for parameter transform purposes
dist.min_vals = min_vals;
dist.max_vals = max_vals;

end

