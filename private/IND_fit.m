function [GMModel, betamix_params] = IND_fit( part_thetas, prior )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   num_components = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read out the number of parameters
[Nparts, Ntheta] = size(part_thetas);
% Initialise the parameters of the beta mixtures used for marginals
betamix_params = zeros( Ntheta, 5 );
% Initialise other properties for transformed particles
U = zeros(Nparts, Ntheta);
Z = zeros(Nparts, Ntheta);
log_marg = zeros(Nparts, 1);
log_Z = zeros(Nparts, 1);

% Transform particle locations for each parameter to [0,1] in order to fit beta distributions
part_thetas_trans = ( part_thetas - prior.theta_min ) ./ ( prior.theta_max - prior.theta_min );

% Fit beta mixtures to each parameter individually, and store
% weight corrections for calculation of logpdf too
for k = 1:Ntheta
    
    % Find a best-fit beta distribution for this parameter
    beta_params = betafit( part_thetas_trans(:,k) );
    
    % Use the best-fit beta as initialisation for a two-component beta mixture
    negbetamix_loglike_here = @(params) -betamix_loglike( part_thetas_trans(:,k), params );
    [temp_params,~] = fminsearch( negbetamix_loglike_here, [0, log(beta_params), log(beta_params)], optimset('Display','off') );
    % Undo the transformation used to allow unconstrained optimisation
    betamix_params(k,1) = 1 / ( 1 + exp(-temp_params(1)) );
    betamix_params(k,2:5) = exp( temp_params(2:5) );
    
    % Find the CDF value of the beta mixture corresponding to each particle
    % These will be approximately uniformly distributed thanks to the probability integral transform
    U(:,k) = betamixcdf( part_thetas_trans(:,k), betamix_params(k,:) );
    
    % Transform this 'uniform' data to be normally distributed
    Z(:,k) = norminv( U(:,k) );
    
    % Calculate the contributions for this parameter to the betamix PDF across all parameters
    log_marg = log_marg + log( betamixpdf( part_thetas_trans(:,k), betamix_params(k,:) ) );
    
    % Calculate the contributions for this parameter to the normal PDF across all parameters
    log_Z = log_Z + log( normpdf( Z(:,k), 0, 1 ) );
    
end

% Fit a Gaussian mixture model to the normal-transformed particle locations
try
    GMModel = fitgmdist(Z, num_components, 'RegularizationValue', 0.01, 'Replicates', 5, 'Options', optimset('Maxiter',1000));
    
% If it fails, just use a one component MVN
catch
    GMModel = fitgmdist(Z, 1, 'Replicates', 5, 'Options', optimset('Maxiter',1000));
    
end