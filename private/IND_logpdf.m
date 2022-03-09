function logpdf = IND_logpdf( theta, GMModel, betamix_params, prior )
% This function calculates the log of the probability density function for
% the copula model created using fitCopulaGMM, with parameters as input

% Read out number of parameters
Nparam = length(prior.theta_min);

% Scale this parameter back to [0,1] using the prior
theta = ( theta - prior.theta_min ) ./ ( prior.theta_max - prior.theta_min );

% Find CDF of the beta mixture for each parameter to get uniform distribution values for each
% Also add the logPDF contributions of each betamix part
log_marg = 0;
U = zeros(1,Nparam);
for k = 1:Nparam
    U(k) = betamixcdf( theta(k), betamix_params(k,:) );
    log_marg = log_marg + log( betamixpdf( theta(k), betamix_params(k,:) ) );
end

% Transform each to be normally distributed
Z = norminv(U);

% Calculate the contribution to logPDF from each component from the normal
log_Z = sum( log( normpdf( Z, 0, 1 ) ) );

% Find the probability that the GMM gave this output Z and add the
log_GMM = log( pdf(GMModel, Z) );

% Combine all factors to get the full logPDF
logpdf = log_GMM - log_Z + log_marg;

end