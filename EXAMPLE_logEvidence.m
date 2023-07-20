function EXAMPLE_logEvidence
% This function demonstrates and validates the use of the SMC routine as a
% means of measuring the log-evidence (logarithm of the posterior
% normalising constant)


%%% Create a test problem consisting of a product of two random Gaussians

% Create first random Gaussian
mu0 = 5 * randn(1,3);
L = 2 * [randn, 0, 0; randn, randn, 0; randn, randn, randn];
SIGMA0 = L * L';
% Create second random Gaussian
mu1 = 5 * randn(1,3);
L = 2 * [randn, 0, 0; randn, randn, 0; randn, randn, randn];
SIGMA1 = L * L';

% Form a posterior that is the product of two Gaussians by assigning one as
% prior, the other as likelihood
prior = distNormal( mu0, SIGMA0 );
like_dist = distNormal( mu1, SIGMA1 );

% Use a model y = theta so that the likelihood just specifies a
% distribution over theta directly
f_model = @(theta) theta;

% Create a "likelihood" function that's just the second Gaussian pdf
f_loglikelihood = @(y,dummy) like_dist.logpdf(y);



%%% Run the SMC

% Define options
options = struct('jumpType', 'MVN', 'Nparts', 10000, 'alpha', 0.75);
% Call performSMC function
[particles, logevidence, options, diagnostics] = performSMC( f_model, prior, f_loglikelihood, options );

% Save the results
save(['OUTPUT_evidencetest','.mat'], 'particles', 'logevidence', 'prior', 'f_model', 'f_loglikelihood', 'options', 'diagnostics');


%%% Calculate the closed form value of the log-evidence and output results

% Pre-calculate important quantities
d = length(mu0);
invSIGMA0 = SIGMA0 \ eye(d);
invSIGMA1 = SIGMA1 \ eye(d);
SIGMA01 = ( invSIGMA0 + invSIGMA1 ) \ eye(d);
mu01 = ( SIGMA01 * ( invSIGMA0 * mu0' + invSIGMA1 * mu1' ) )';

% Closed form of log-evidence, obtained by comparing the raw product of two
% Gaussian pdfs with the correctly-normalised Gaussian pdf
log_Z = -0.5 * d * log(2*pi) - 0.5 * log( det(SIGMA0 * SIGMA1) ) + 0.5 * log( det(SIGMA01) ) + 0.5 * ( mu01 * (SIGMA01 \ mu01') - mu0 * (invSIGMA0 * mu0') - mu1 * (invSIGMA1 * mu1') );

% Print the comparison for the user
fprintf('Actual log-evidence is %g, while SMC-generated estimate is %g\n',log_Z,logevidence);

end