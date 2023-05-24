function w = weightParticles( loglikes, dT )
% This function calculates normalised weights for all of the particles with
% the provided log likelihoods, assuming an increase in temperature of dT

% Weight particles according to loglikelihood
w = dT * loglikes;
   
% Shift all weights to be near zero
w = w - max(w);

% Change to weights based on likelihood (instead of loglikelihood)
w = exp(w);

% Normalise the weights
w = w / sum(w);   

end