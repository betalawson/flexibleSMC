function y = betapdf_raw(x,a,b)
% A version of MATLAB's BETAPDF that has no input handling or error
% checking, for speed reasons.

% Compute logs
smallx = x<0.1;

loga = (a-1).*log(x);
logb = zeros(size(x));
logb(smallx) = (b-1) .* log1p(-x(smallx));
logb(~smallx) = (b-1) .* log(1-x(~smallx));

y = exp(loga+logb - betaln(a,b));
