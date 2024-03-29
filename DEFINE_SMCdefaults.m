function options = DEFINE_SMCdefaults
% This function is where the user can set the options for the SMC

options.Nparts = 5000;                  % Number of particles
options.alpha = 0.5;                    % Keep ratio
options.minMCMCsteps = 5;               % Minimum number of MCMC steps (number to perform before estimating number of steps to take)
options.maxMCMCsteps = 100;             % Maximum number of MCMC steps

options.jumpType = 'MVN';               % Type of jumping distribution to use
options.meanFlip = 2;                   % Maximum number of bits to flip for binary updates

%%%%% SMC-ABC specific %%%%%
options.discrepancy = 'mahalanobis';    % Discrepancy function - may provide a function @(S) that returns the discrepancy for a set of summary statistics, or may specify 'euclid', 'weighted', 'mahalanobis', 'SSnonans'
options.mutateAll = 'true';             % Specifies whether all particles will have their positions updated by the mutation step (used in SMC-ABC)
options.D_target = -Inf;                % Target discrepancy that triggers exit from the SMC-ABC routine. Set to -Inf to continue reducing discrepancy threshold until particle uniqueness lost

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.verbose = true;                     % Display output while running
options.visualise = true;                   % Display graphics while running
options.visFunc = @defaultVisualisation;    % Function to use for graphical display

end