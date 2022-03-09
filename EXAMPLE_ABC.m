function EXAMPLE_ABC


%%% Initial specifications

% Specify prior
prior = priorUniform( [0.5*ones(1,11), 0.1], [1.5*ones(1,11), 50] );
% Specify model
f_model = @(theta) modelCRN(theta);


%%% Generate some synthetic data

% Specify the true parameter values
theta_true = [ones(1,11)];
% Specify the noise level
sigma = 2;

% Merge together to be a single parameter
theta_true = [theta_true, sigma];

% Create noisy observations from the model
target_data = f_model(theta_true);
target_data = target_data + sigma .* randn(size(target_data));


%%% Specify likelihood (may depend on data)
f_loglikelihood = @(y, theta) likelihoodGaussianVarySigma(y, target_data, theta);


%%% Run the SMC

% Define options
options = struct('jumpType', 'IND', 'Nparts', 500);
% Call performSMC function
[particles, diagnostics] = performSMC( f_model, prior, f_loglikelihood, options );

% Save the results
save(['test_output','.mat'], 'particles', 'diagnostics');

end