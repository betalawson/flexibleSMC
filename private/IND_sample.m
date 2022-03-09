function theta = IND_sample( GMModel, betamix_params, prior )
% This function samples from a copula model in which parameter marginals
% are transformed to be approximately normal after fitting a mixture of two
% betas to their original [0,1] transformed marginals. The transformed
% parameters are fit together using a mixture of Gaussians.

% Read out number of parameters
Ntheta = length(prior.theta_min);

% Pick a new location from the Gaussian mixture model
Z_prop = random( GMModel );

% Convert this to a uniformly distributed random variable (in each parameter)
U_prop = normcdf( Z_prop, 0, 1 );

% Loop over each parameter, converting it from its uniform value back to
% its original 
theta = zeros(1,Ntheta);
for j = 1:Ntheta
    betamixqantfun = @(u) betamixcdf(u, betamix_params(j,:) ) - U_prop(j);
    theta(j) = bisectionSolve( betamixqantfun, [0,1] );
end

% Undo the scaling to [0, 1]
theta = prior.theta_min + (prior.theta_max - prior.theta_min) .* theta;