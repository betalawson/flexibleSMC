function theta = REG_sample(X,Y,Y_obs,T,C_params,C_sigmas,temp)

% Read out dimensions
D_params = size(C_params,2);
D_sigmas = size(C_sigmas,2);

% Sample in the parameter space
params = X(1:D_params) + temp*(Y_obs(:) - Y(:))' * T + mvnrnd(zeros(D_params,1), C_params);

%    Current +  Linear Correction   +    MVN including nullspace adjustment

% Sample in the sigma parameter space
sigma_params = X(D_params+1:end) + mvnrnd(zeros(D_sigmas,1), C_sigmas);

% Combine together to get a parameter vector
theta = [params, sigma_params];

end

