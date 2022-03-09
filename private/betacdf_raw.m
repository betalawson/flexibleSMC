function p = betacdf_raw(x,a,b)
% This is an internal function based on MATLAB's BETACDF function, but with
% no smart handling of inputs or checking of the distribution for validity.
% Use of this function matches the basic use of BETACDF(x,a,b).

% Weed out any out of range parameters or data.
okAB = (0 < a & a < Inf) & (0 < b & b < Inf);
k = (okAB & (0 <= x & x <= 1));
allOK = all(k(:));

% Fill in NaNs for out of range cases, fill in edges cases when X is outside 0 or 1.
if ~allOK
    p = NaN(size(k),'like',internal.stats.dominantType(x,a,b)); % Single if x, a, or b is
    p(okAB & x < 0) = 0;
    p(okAB & x > 1) = 1;
    
    % Remove the bad/edge cases, leaving the easy cases.  If there's
    % nothing remaining, return.
    if any(k(:))
        if numel(x) > 1, x = x(k); end
        if numel(a) > 1, a = a(k); end
        if numel(b) > 1, b = b(k); end
    else
        return;
    end
end

pk = betainc(x,a,b);

% Broadcast the values to the correct place if need be.
if allOK
    p = pk;
else
    p(k) = pk;
end
