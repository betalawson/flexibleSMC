function pdf = betamixpdf( x, betamix_params )
% Evaluates the probability density function for the two-component beta
% mixture with parameters defined by (betamix_params) at the point
% indicated by (x)

pdf = betamix_params(1) * betapdf_raw( x, betamix_params(2), betamix_params(3) ) + (1 - betamix_params(1)) * betapdf_raw( x, betamix_params(4), betamix_params(5) );

end

