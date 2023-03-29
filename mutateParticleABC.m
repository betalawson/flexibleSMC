function [P, accepted, move_runs] = mutateParticleABC( P, J, epsilon, prior, f_model, f_summaries, f_discrep )
% This function mutates particles according to an ABC "likelihood". If a
% particle's discrepancy value falls within the current allowable threshold
% 'epsilon', it has a "likelihood" of one, otherwise zero.

% Perform each move sequentially
accepted = false;
move_runs = 0;
    
% Propose new location according to jumping distribution J
prop_theta = J.sample( P.theta );
    
% Check prior in this location
prop_logprior = prior.logpdf( prop_theta );
    
% Evaluate particle if it's within the prior
if prop_logprior > -Inf
    
    % Calculate the Metropolis-Hastings ratio
    if J.symmetric
        MH = exp( prop_logprior - P.logprior );
    else
        MH = exp( prop_logprior - P.logprior + J.logpdf( P.theta, prop_theta ) - J.logpdf( prop_theta, P.theta ) );
    end
    
    % Check if the particle has been rejected by the M-H ratio before doing any simulation
    if rand < min( [1,MH] )
        
        % Run the model and calculate associated summaries and discrepancy
        prop_y = f_model( prop_theta );
        prop_S = f_summaries( prop_y );
        prop_D = f_discrep( prop_S );
        move_runs = move_runs + 1;
        
        % Accept this move if discrepancy falls under threshold
        if prop_D <= epsilon
            P = struct('theta', prop_theta, 'y', prop_y, 'logprior', prop_logprior, 'S', prop_S, 'D', prop_D);
            accepted = true;
        end
        
    end
end