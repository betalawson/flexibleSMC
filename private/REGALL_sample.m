function theta = REGALL_sample(X,Y,Y_obs,T,C,temp)

% Read out dimensions
theta = X + temp * (Y_obs(:) - Y(:))' * T + mvnrnd(zeros(1,size(C,2)), C);


end

