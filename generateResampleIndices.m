function r = generateResampleIndices(w, method)
% This function takes as input a vector of N weights, and generates N
% samples (indices associated with the weight vector) that satisfy (in
% expectation) that each index is sampled a number of times according to
% its weight,
%    E[ sum_i I(r_i = j) ] = w_j N.
%
% Resampling method 'multinomial' is a common way to achieve this, but the
% SMC literature also recommends other approaches that achieve less
% variance in the number of indices assigned to each particle:
%
%    'multinomial':  generate N assignments independently at random
%    'stratified':   generate 1 assignment randomly from each bin of width
%                    1/N along the "cdf" of weights
%    'residual':     generate N-M assignments deterministically, that match
%                    the rounded-down values w_j N, then use 'stratified'
%                    for the remaining M
%    'residual-multinomial':    as in 'residual', but the M remaining
%                               assignments use 'multinomial'

% Count the number of particles given using length of weight vector
N = length(w);

% Method parsing
switch lower(method)
    
    case 'residual'
        deterministic = true;
        remainder_method = 'stratified';
        
    case 'residual-multinomial'
        deterministic = true;
        remainder_method = 'multinomial';
        
    case 'stratified'
        deterministic = false;
        remainder_method = 'stratified';
        
    case 'multinomial'
        deterministic = false;
        remainder_method = 'multinomial';
        
    otherwise
        error('The requested resampling method has not been implemented (or entered incorrectly!)');
        
end

% Immediately normalise the weights
w = w / sum(w);

% If a 'residual' based method is being used, perform the deterministic
% component of the assignment first
if deterministic
    
    % Determine guaranteed counts and assign indices
    expected_counts = w*N;
    guarantee_counts = floor(expected_counts);
    r_fixed = repelem(1:N, guarantee_counts);
    % Count how many assignments are left to perform
    N_random = N - length(r_fixed);
    % Use the "leftover" weighting as the new weights and re-normalise
    w = expected_counts - guarantee_counts;
    w = w / sum(w);
    
else
    
    % No fixed assignments, all assignments are random
    r_fixed = [];
    N_random = N;
    
    
end

% Now, perform N_random random assignments using the remainder_method
if N_random > 0

    switch remainder_method
    
        case 'multinomial'
            % Each individual assignment is random sampling with replacement
            r_random = randsample(1:N, N_random, true, w);
        
        case 'stratified'
               
            % Generate random positions along the CDF to sample from
            get_locs = ( (0:N_random-1) + rand(1,N_random) ) / N_random;
            % Get the cdf for weights, ensure last value is above one
            w_cdf = cumsum(w);
            w_cdf(end) = 1 + sqrt(eps);
            
            % Now grab out the "get_locs" locations 
            r_random = discretize(get_locs, [0 w_cdf]);
            
    end
    
else
   
    % No random assignments performed, so set it to blank
    r_random = [];
    
end

% Concatenate together the random and fixed assignments
r = [r_fixed, r_random];

end