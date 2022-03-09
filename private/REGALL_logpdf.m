function logpdf = REGALL_logpdf(X,X_old,Y,Y_obs,T,C,temp)

logpdf =   mvnlogpdf( X, X_old + temp*(Y_obs(:) - Y(:))' * T, C );


     
end

