function [X, W] = whitenData(X)
% This function takes an input data matrix X, and returns it as a whitened
% data matrix, along with the transform matrix used to do so

% Subtract the mean of the data
X = X - mean(X);

% Find the eigendecomposition of the covariance matrix, adding minor
% diagonal stabilisation
stabiliser = -17;
found = false;
while ~found
   
    W = ( cov(X) + 10^(stabiliser) * eye(size(X,2)) )^(-1/2);
    if all(isreal(W),'all')
        found = true;
    else
        stabiliser = stabiliser + 1;
        if stabiliser > -10
            warning('Using stabilisation larger than 1e-10');
        end
    end
    
end

% Update eigenvalues and eigenvectors based on this better behaved matrix
%W = LAM^(-1/2) * V';

% Apply it to the data
X = (W * X')';

end

