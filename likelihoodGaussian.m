function loglike = likelihoodGaussian( y, y_obs, sigma )
% This function calculates a plain Gaussian (log) likelihood between the
% model output y, and observations y_obs. The standard deviation may be
% specified as a vector/matrix of appropriate dimensions to match y, or as 
% a single value.

% Log likelihood is found by summing the term contributed by each
% individual datapoint
if any( isnan(y) )
    loglike = -Inf;
else
    loglike = sum( -log(sigma) - 0.5 * log( 2 * pi ) - 0.5 * (y - y_obs).^2 ./ sigma.^2, 'all' );
end

end

