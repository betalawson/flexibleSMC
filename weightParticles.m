function w = weightParticles( loglikes, dT )
% This function calculates weights for all of the particles with the
% provided log likelihoods, assuming an increase in temperature of dT

% Weight particles according to loglikelihood
w = dT * loglikes;

% Shift all weights to be near zero
w = w - max(w);

% Weights based on actual likelihood
w = exp(w);

end