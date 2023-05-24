function EXAMPLE_logEvidence

%%% Use a product of random Gaussians to get a closed form for the log-evidence

% Define random mean vectors and covariance matrices
mu0 = 5 * randn(3,1);
mu1 = 5 * randn(3,1);
L = 2 * [randn, 0, 0; randn, randn, 0; randn, randn, randn];
SIGMA0 = L * L';
L = 2 * [randn, 0, 0; randn, randn, 0; randn, randn, randn];
SIGMA1 = L * L';

% Specify prior
prior = distNormal( mu0, SIGMA0 );
% Specify "likelihood" distribution
like_dist = distNormal( mu1, SIGMA1 );

% Create a dummy model that just returns the same value
f_model = @(theta) theta;

% Create a "likelihood" function that's just the second Gaussian pdf
f_loglikelihood = @(y,dummy) like_dist.logpdf(y);



%%% Run the SMC

% Define options
options = struct('jumpType', 'MVN', 'Nparts', 100000, 'alpha', 0.75);
% Call performSMC function
[particles, logevidence, options, diagnostics] = performSMC( f_model, prior, f_loglikelihood, options );

% Save the results
save(['test_output','.mat'], 'particles', 'logevidence', 'prior', 'f_model', 'f_loglikelihood', 'options', 'diagnostics');



%%% Compare the logevidence estimate to the closed-form value
log_Z = 0.5 * length(mu0) * log(2*pi) + 0.5 * log( det(SIGMA0 + SIGMA1) ) + 0.5 * (mu1 - mu0)' * ( (SIGMA0 + SIGMA1) \ (mu1 - mu0) ) + 0.5 * length(mu0) * log(2*pi) - 0.5 * log(det( (SIGMA0\eye(length(mu0))) + (SIGMA1\eye(length(mu1))) ));
fprintf('Actual ''''log-evidence'''' is %g, while SMC-generated estimate is %g\n',log_Z,logevidence);

end