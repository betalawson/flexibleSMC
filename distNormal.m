function dist = distNormal(mu, SIGMA)

% Create a sampling function
dist.sample = @() mvnrnd(mu, SIGMA);

% Create a log pdf evaluation function
U = cholcov(SIGMA);
const = -0.5 * length(mu) * log(2*pi) - sum( log(diag(U)) );
dist.logpdf = @(x) const - 0.5 * norm( U' \ (x - mu)', 2)^2;