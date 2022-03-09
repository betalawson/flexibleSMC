function [particles, diagnostics] = performSMC( f_model, prior, f_loglikelihood, target_data, input_options )
% This function performs SMC-based inference on a supplied model, with
% user-specified functions for the prior and the likelihood. Matching to
% specific target data should be baked into the likelihood function - see
% the attached example for this procedure in action.



%%% Option initialisation

% Load default options and adjust as user has requested
options = DEFINE_SMCdefaults();
provided_options = fieldnames(input_options);
for k = 1:length(provided_options)
    options.(provided_options{k}) = input_options.(provided_options{k});
end

% Read out number of particles
Nparts = options.Nparts;



%%% Create the initial particles
parfor k = 1:Nparts
   
    % Loop until a suitable particle is found
    particle_found = false;
    prior_runs(k) = 0;
    while ~particle_found
    
        % Sample from the prior
        prop_theta = prior.sample();
    
        % Simulate the model here
        prop_y = f_model( prop_theta );
        prior_runs(k) = prior_runs(k) + 1;
    
        % Evaluate the likelihood of this particle 
        prop_loglike = f_loglikelihood( prop_y, prop_theta ); 
    
        % Only accept profiles with nonzero likelihood
        if prop_loglike > -Inf
            particle_found = true;
            particles{k} = struct('theta', prop_theta, 'y', prop_y, 'logprior', prior.logpdf( prop_theta ), 'loglike', prop_loglike);
        end
        
    end
    
end

% Store total runs used to generate the particles, output if verbose
model_runs = sum(prior_runs);
if options.verbose
    fprintf('Prior was successfully sampled, with %g%% of prior samples accepted\n', Nparts / model_runs * 100 );
end



%%% Basic Setup

% Create a function that mutates particles in terms of only the varying properties
mutateParticle = @(P, J, T) mutateParticleFull( P, J, T, prior, f_model, f_loglikelihood );

% Initialise acceptance rate vector
acceptance_rates = [];



%%% Visualise if Requested
if options.visualise
    options.visFunc( particles );
    pause(0.5);
end


%%% Perform the resample-move SMC algorithm

% Initialise temperature at zero
T = 0;

% Increment temperature until the posterior (T = 1) is being sampled from
while T < 1
    
    
    %%% SELECT NEW TEMPERATURE
    
    % Grab log likelihoods of all particles to perform temp. calculations
    part_loglikes = cellfun( @(x) x.loglike, particles );
    
    % First trial a temperature of one to see if its ESS is satisfactory
    trial_ESS = calculateESS( weightParticles( part_loglikes, 1-T ) );
    
    % Accept this T value if acceptable, or use bisection to find highest
    % possible T that does meet the desired ESS threshold
    if trial_ESS >= options.alpha * Nparts
        T_new = 1;
    else
        ESS_func = @(T_new) calculateESS( weightParticles( part_loglikes, T_new - T ) ) - options.alpha * Nparts;
        T_new = bisectionSolve( ESS_func, [T+eps 1] );
    end
    
    % Calculate particle weights using the temperature found
    part_ws = weightParticles( part_loglikes, T_new - T );
    
    
    %%% RESAMPLE PARTICLES
    
    % Weighted resampling using particles' loglike-derived weights
    r = randsample(1:Nparts, Nparts, true, part_ws);
    % Copy particles using the weighted resampling
    particles = particles(r);
    
    
    %%% SET UP JUMPING DISTRIBUTION
    
    % Use the constructor function for the jumping distribution
    J = constructJumpingDistribution( particles, prior, target_data, T_new, options );
    
    
    
    %%% MUTATE (MOVE) PARTICLES
    
    % Grab out the minimum number of jumps just for code cleanliness
    R_min = options.minMCMCsteps;
    
    % Perform an initial number of MCMC steps to estimate acceptance rate,
    % used to calculate how many MCMC steps will be taken     
    accept_rates = zeros(1,R_min);
    for k = 1:R_min
        for m = 1:Nparts
            [ particles{m}, accepted(m), move_runs(m) ] = mutateParticle( particles{m}, J, T_new );
        end
        accept_rates(k) = mean(accepted);
    end
    model_runs = model_runs + sum( move_runs );
    
    % Estimated acceptance rate is the mean across all particles
    est_acc_rate = mean( accept_rates );
    % Force acceptance rate into the range [1e6, 1 - 1e6]
    est_acc_rate_trim = max( [ min( [est_acc_rate, 1-1e-6] ), 1e-6 ] );
    % Calculate number of additional steps to perform
    R = ceil( log(0.05) / log( 1 - est_acc_rate_trim ) );
    % Ensure recommended number of steps does not exceed maximum
    R = min( [R, options.maxMCMCsteps] );
    % Convert this into a number of additional runs to perform
    R_add = max( [R - R_min, 0] );
    
    % Perform the remaining steps as suggested by recommended stepcount R
    for k = 1:R_add
        for m = 1:Nparts
            [ particles{m}, accepted(m), move_runs(m) ] = mutateParticle( particles{m}, J, T_new );
        end
        accept_rates(k) = mean(accepted);
    end
    model_runs = model_runs + sum( move_runs );  
    
    
    %%% UPDATE TEMPERATURE AND DIAGNOSTICS
    
    T = T_new;
    overall_acceptance = ( R_min * est_acc_rate + R_add * mean( accept_rates ) ) / ( R_min + R_add );
    acceptance_rates(end+1:end+R_min+R_add) = overall_acceptance;
    if options.verbose
        fprintf('Step taken using T = %g, with acceptance rate %g%%\n', T, overall_acceptance * 100);
    end
    
    
    %%% VISUALISE IF REQUESTED
    if options.visualise
        options.visFunc( particles );
        pause(0.5);
    end
    
end



%%% FINAL STORAGE

% Save diagnostic information
diagnostics.model_runs = model_runs;
diagnostics.acceptance_rates = acceptance_rates;
diagnostics.overall_acceptance_rate = mean(acceptance_rates);       