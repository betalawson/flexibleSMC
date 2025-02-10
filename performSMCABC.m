function [particles, options, diagnostics] = performSMCABC( f_model, f_summaries, S_target, prior, input_options )
% This function performs SMC-ABC-based inference on a supplied model, with
% the user also providing a function that generates summary statistics.
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
for k = 1:Nparts
   
    % Loop until a suitable particle is found
    particle_found = false;
    prior_runs(k) = 0;
    while ~particle_found
    
        % Sample from the prior
        prop_theta = prior.sample();
    
        % Simulate the model here
        prop_y = f_model( prop_theta );
        prior_runs(k) = prior_runs(k) + 1;
        
        % Summarise the model output
        prop_S = f_summaries( prop_y );
       
        % If particle is acceptable, keep it
        if ~any( isnan(prop_S) )
            particle_found = true;
            particles{k} = struct('theta', prop_theta, 'y', prop_y, 'logprior', prior.logpdf( prop_theta ), 'S', prop_S);
        end
        
    end
    
end

% Store total runs used to generate the particles, output if verbose
model_runs = sum(prior_runs);
if options.verbose
    fprintf('Prior was successfully sampled, with %g%% of prior samples accepted\n', Nparts / model_runs * 100 );
end


%%% Setup the discrepancy function

% If discrepancy is given as a string, create the function for discrepancy
if isa(options.discrepancy,'char')
    switch lower(options.discrepancy)
        
        % Normal Euclidean distance
        case {'euclid','euclidean','standard'}
            f_discrep = @(S) pdist2(S, S_target, 'euclidean');
            
        % Weighted Euclidean using particle variances
        case {'weighted', 'weight', 'weightedeuclidean','weighted-euclidean'}
            f_discrep = @(S) pdsit2(S, S_target, 'seuclidean', 'DistParameter', sqrt( diag( SIGMA0 ) ) );
        
        % Weighted Euclidean that removes summary correlations
        case {'mahalanobis'}
            f_discrep = @(S) pdist2(S, S_target, 'mahalanobis', 'DistParameter', SIGMA0 );
            
        % Sum of squares, ignoring NaNs
        case {'ssnonans'}
            f_discrep = @(S) sum( (S - S_target).^2, 2, 'omitnan');

    end
    
% If discrepancy is given as a function, use it to define f_discrep    
elseif isa(options.discrepancy,'function_handle')
    f_discrep = @(S) options.discrepancy(S, S_target);
    
else
    error('Discrepancy function not specified correctly!');
end


%%% Other basic Setup

% Create a function that mutates particles in terms of only the varying properties
mutateParticle = @(P, J, epsilon) mutateParticleABC( P, J, epsilon, prior, f_model, f_summaries, f_discrep );

% Determine the discrepancies of all particles
part_Ds = f_discrep( getProperty(particles,'S') );
% Set these values in the particle struct
particles = addProperty( particles, 'D', part_Ds );

% Determine the number of the worst particle to keep
Nkeep = floor( Nparts * options.alpha );
% According to user options, create the list of which particles to mutate
if options.mutateAll
    mutate_range = 1:Nparts;
else
    mutate_range = Nkeep+1:Nparts;
end

% Initialise diagnostic storage
acceptance_rates = [];



%%% Perform the resample-move style SMC-ABC algorithm

looping = true;
while looping
    
    
    %%% SORT PARTICLES BY DISCREPANCY
    
    % Get particle discrepancies
    part_Ds = getProperty(particles, 'D');
    
    % Sort particles using these
    [~,I] = sort(part_Ds, 'ascend');
    particles = particles(I);
    
    
    
    %%% DETERMINE CULL THRESHOLD AND RESAMPLE REJECTED PARTICLES
    
    % Find threshold by checking worst kept particle
    D_threshold = particles{Nkeep}.D;
    
    % Choose a random sample from kept particles to copy onto
    r = randi(Nkeep, [Nparts - Nkeep, 1]);
    particles(Nkeep+1:Nparts) = particles(r);
    
    
    
    %%% SET UP JUMPING DISTRIBUTION
    
    % Use the constructor function for the jumping distribution
    J = constructJumpingDistribution( particles, prior, options );
    
    
    
    %%% MUTATE (MOVE) PARTICLES
    
    % Grab out the minimum number of jumps just for code cleanliness
    R_min = options.minMCMCsteps;
    
    % Perform an initial number of MCMC steps to estimate acceptance rate,
    % used to calculate how many MCMC steps will be taken     
    accept_rates = zeros(1,R_min);
    for k = 1:R_min
        accepted = zeros(1,Nparts);
        move_runs = zeros(1,Nparts);
        parfor m = mutate_range
            [ particles{m}, accepted(m), move_runs(m) ] = mutateParticle( particles{m}, J, D_threshold );
        end
        accept_rates(k) = mean(accepted(mutate_range));
        model_runs = model_runs + sum( move_runs );
    end
    
    % Estimated acceptance rate is the mean across all particles
    est_acc_rate = mean( accept_rates );
    % Force acceptance rate into the range [1e6, 1 - 1e6]
    est_acc_rate_trim = max( [ min( [est_acc_rate, 1-1e-6] ), 1e-6 ] );
    % Calculate number of additional steps to perform
    R_rec = ceil( log(0.05) / log( 1 - est_acc_rate_trim ) );
    % Ensure recommended number of steps does not exceed maximum
    R = min( [R_rec, options.maxMCMCsteps] );
    % Convert this into a number of additional runs to perform
    R_add = max( [R - R_min, 0] );
    
    % Perform the remaining steps as suggested by recommended stepcount R
    for k = 1:R_add
        move_runs = zeros(1,Nparts);
        parfor m = mutate_range
            [ particles{m}, accepted(m), move_runs(m) ] = mutateParticle( particles{m}, J, D_threshold );
        end
        accept_rates(k) = mean(accepted);
        model_runs = model_runs + sum( move_runs );  
    end
    
    
    %%% UPDATE DIAGNOSTICS
    
    overall_acceptance = ( R_min * est_acc_rate + R_add * mean( accept_rates ) ) / ( R_min + R_add );
    acceptance_rates(end+1:end+R_min+R_add) = overall_acceptance;
    if options.verbose
        fprintf('Step taken, acceptable discrepancy was D = %g. Acceptance rate was %g%%\n', D_threshold, overall_acceptance * 100);
    end
    
    
    %%% LOOP TERMINATION
    
    % Loop terminates if the target discrepancy is reached
    if D_threshold <= options.D_target
        fprintf('SMC-ABC routine finished: target discrepancy reached.\n');
        looping = false;
    end
    
    % Loop terminates if acceptance rate fell too low and R exceeded max. R
    if R_rec > options.maxMCMCsteps
        fprintf('SMC-ABC routine terminated: acceptance rate too low.\n');
        looping = false;
    end
    
end


%%% FINAL STORAGE

% Save diagnostic information
diagnostics.model_runs = model_runs;
diagnostics.acceptance_rates = acceptance_rates;
diagnostics.overall_acceptance_rate = mean(acceptance_rates);       