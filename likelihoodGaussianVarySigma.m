function loglike = likelihoodGaussianVarySigma( y, y_obs, theta )
% This function calculates a plain Gaussian (log) likelihood between the
% model output y, and observations y_obs. The standard deviation may be
% specified as a vector the same length as y, or as a single value.

% Final elements of theta are the standard deviations of error term
% associated with each individual component of y.
Ny = size(y,1);
sigma = theta(end-Ny+1:end);

% Log likelihood is found by summing the term contributed by each
% individual datapoint
if any( isnan(y) )
    loglike = -Inf;
else
    loglike = sum( -log(sigma) - 0.5 * log( 2 * pi ) - 0.5 * (y - y_obs).^2 ./ sigma.^2, 'all' );
end

end

