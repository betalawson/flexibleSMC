function J = constructJumpingDistribution( particles, prior, options )
% This function constructs a jumping distribution of the requested type.
% The current particle locations, and the prior, are used.


% Grab out a matrix of all particle theta values
part_thetas = getProperty( particles, 'theta' );

% Construct requested distribution
switch lower(options.jumpType)
       
        
        % Multivariate random normal jumps
        case {'mvn','normal'}
            
            % Prepare jump covariance
            C = 2.38^2 / size(part_thetas,2) * cov(part_thetas);
            % Multivariate normal random jumps, informed by particle covariance
            J.sample = @(theta, ~) mvnrnd( theta, C );
            % Jumping distribution is symmetric so logpdf not needed
            J.symmetric = true;
            
            
        % Independent proposals using a copula approach     
        case {'independent','ind'}
        
            % Find the copula and best-fit marginals for current particle locations
            [ GMModel, betamix_params ] = IND_fit( part_thetas, prior );

            % Use the functions for sampling and logpdf for these 
            J.sample = @(theta, ~) IND_sample( GMModel, betamix_params, prior );
            J.logpdf = @(theta, theta_old, ~) IND_logpdf( theta, GMModel, betamix_params, prior );

            % Mark this as an assymetric jumping distribution because
            %         J(theta'|theta) = J(theta')
            J.symmetric = false;            
        
end