function y_model = modelCRN(params)

% Run the model and get its output
[t_model, y_model] = simulateCourtemanche(params);

% Calculate the action potential biomarkers as a way to check that a
% reasonable AP was obtained
biomarkers = computeBiomarkers(t_model, y_model);

% Reshape model output to be compatible with the remaining code
y_model = y_model';

% If the model solution failed, then make the results ridiculous
if any(isnan(biomarkers))
    y_model = nan(size(y_model));
end

