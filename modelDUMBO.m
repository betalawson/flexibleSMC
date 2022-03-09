function y_model = modelDUMBO(params)

% Just a random linear model to test things
y_model = 100 + params(1:6) * [ [2;3;1;-4;1;-4], [1;-4;2;1;0;2], [0;1;5;2;3;7], [-1;2;-1;2;2;-1], [0;0;1;0;1;0]];