function logpdf = REG_logpdf(X,X_old,Y,Y_obs,T,C_params,C_sigmas,temp)

% Read out dimensions
D_params = size(C_params,2);

% Probability is the product of the probability of generating this
% parameter location, and these values of the sigma parameters. This
% corresponds to addition in terms of logarithms
logpdf =   mvnlogpdf( X(1:D_params), X_old(1:D_params) + temp*(Y_obs(:) - Y(:))' * T, C_params ) ...
         + mvnlogpdf( X(D_params+1:end), X_old(D_params+1:end), C_sigmas );


     
end

