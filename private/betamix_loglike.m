function loglike = betamix_loglike(x, trans_betamix_params)
% This function evaluates the likelihood of observing the data x given a
% beta mixture with parameters betamix_params. However, parameters are
% given in a transformed form so as to reflect the conditions on their
% values

% De-transform the inputs
betamix_params(1) = 1 / ( 1 + exp(- trans_betamix_params(1) ) );
betamix_params(2:5) = exp( trans_betamix_params(2:5) );

% Calculate log likelihood of all data together (sum of logs is log of their product)
loglike = sum( log( betamixpdf(x, betamix_params) ) );