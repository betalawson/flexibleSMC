function [P, accepted, move_runs] = mutateParticleFull( P, J, T, prior, f_model, f_loglikelihood )

% Perform each move sequentially
accepted = false;
move_runs = 0;
    
% Propose new location according to jumping distribution J
prop_theta = J.sample( P.theta, P.y );
    
% Check prior in this location
prop_logprior = prior.logpdf( prop_theta );
    
% Evaluate particle if it's within the prior
if prop_logprior > -Inf
        
    % Evaluate model and log likelihood
    prop_y = f_model( prop_theta );
    prop_loglike = f_loglikelihood( prop_y, prop_theta );
    move_runs = move_runs + 1;
    
    % Calculate the Metropolis-Hastings ratio
    if J.symmetric
        MH = exp( prop_logprior - P.logprior + T*(prop_loglike - P.loglike) );
    else
        MH = exp( prop_logprior - P.logprior + T*(prop_loglike - P.loglike) + J.logpdf( P.theta, prop_theta, P.y ) - J.logpdf( prop_theta, P.theta, prop_y ) );
    end
            
    % Accept this particle according to Metropolis-Hastings ratio
    if rand < min( [1,MH] ) && prop_loglike > -Inf
        P = struct('theta', prop_theta, 'y', prop_y, 'logprior', prop_logprior, 'loglike', prop_loglike );
        accepted = true;
    end
    

    
end