function EXAMPLE


%%% Initial specifications

% Specify prior
prior = priorUniform( [0.5*ones(1,6), 0.001], [1.5*ones(1,6), 10] );
% Specify model
f_model = @(theta) modelLinear(theta);


%%% Generate some synthetic data

% Specify the true parameter values
theta_true = [0.6, 1.3, 1.1, 0.8, 0.8, 0.6];
% Specify the noise level
sigma = 0.01;


% Merge together to be a single parameter
theta_true = [theta_true, sigma];

% Create noisy observations from the model
target_data = f_model(theta_true);
target_data = target_data + sigma .* randn(size(target_data));

% Count the number of sigma
Nsigmas = size(target_data,1);

%%% Specify likelihood (may depend on data)
f_loglikelihood = @(y, theta) likelihoodGaussianVarySigma(y, target_data, theta);


%%% Run the SMC

% Define options
options = struct('jumpType', 'MVN', 'Nparts', 1000, 'Nsigmas', Nsigmas, 'alpha', 0.75);
% Call performSMC function
[particles, logevidence, options, diagnostics] = performSMC( f_model, prior, f_loglikelihood, options );

% Save the results
save(['test_output','.mat'], 'particles', 'logevidence', 'prior', 'target_data', 'f_model', 'f_loglikelihood', 'options', 'diagnostics');

end