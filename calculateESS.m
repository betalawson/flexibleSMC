function ESS = calculateESS(w)
% A function that estimates the sample size using the particle weights

ESS = 1 / sum(w.^2);

end

