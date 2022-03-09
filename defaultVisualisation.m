function defaultVisualisation(particles)
% This is a default function that visualises the current particles. Custom
% visualisation functions specific to your problem should probably be
% provided, but this serves as a placeholder and a template for your own
% functions.

%%% READ OUT DATA
thetas = getProperty(particles, 'theta');
Nthetas = size(thetas,2);

%%% DEFINITIONS
theta_names = cell(1,Nthetas);
for k = 1:Nthetas
    theta_names{k} = ['\theta_{',num2str(k),'}'];
end


%%% PARAMETER MARGINALS

% Create figure for parameter marginals if it doesn't exist
figBivariate = findobj('type','Figure','Name','Parameter Marginals');
if isempty(figBivariate)
    figure('units','Normalized','OuterPosition',[0 0 1 1], 'Name', 'Parameter Marginals');
else
    figure(figBivariate);
end

% Clear the figure's current data
clf;

% Figure setup
Ncols = ceil( sqrt(2.5*Nthetas) );     % Approx. twice as many rows as columns
Nrows = ceil( Nthetas / Ncols );     % Rows as necessary to contain plot
ax = createAxes( Nthetas, Nrows, Ncols, [], true );

% Plot parameter marginals
for k = 1:Nthetas
    [f,t] = ksdensity( thetas(:,k) );
    plot(ax{k}, t, f, 'LineWidth', 2);
    xlabel(ax{k}, theta_names{k}, 'FontSize', 20);
    set(ax{k},'FontSize',18);
    set(ax{k},'YTick',[]);
end


%%% PARAMETER BIVARIATE SCATTERS

% Create figure for parameter bivariate plots if it doesn't exist
figBivariate = findobj('type','Figure','Name','Bivariate Scatters');
if isempty(figBivariate)
    figure('units','Normalized','OuterPosition',[0 0 1 1], 'Name', 'Bivariate Scatters');
else
    figure(figBivariate);
end

% Clear the figure's current data
clf;

% Figure setup
Nfigs = Nthetas * (Nthetas - 1) / 2;
Ncols = ceil( sqrt(2.5*Nfigs) );
Nrows = ceil( Nfigs / Ncols );
ax = createAxes( Nfigs, Nrows, Ncols, [], true );

% Loop over all parameter combinations
c = 0;
for i = 1:Nthetas-1
    for j = i+1:Nthetas

        % Increment counter
        c = c + 1;
        % Plot scatter
        plot(ax{c}, thetas(:,i), thetas(:,j), 'k.', 'MarkerSize', 15 );
        xlabel(ax{c}, theta_names{i}, 'FontSize', 16);
        ylabel(ax{c}, theta_names{j}, 'FontSize', 16);
        set(ax{c}, 'FontSize', 16);
    end
end


end