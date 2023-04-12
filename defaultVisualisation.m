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
figObj = findobj('type','Figure','Name','Parameter Marginals');
if isempty(figObj)
    figure('units','Normalized','OuterPosition',[0 0 1 1], 'Name', 'Parameter Marginals');
else
    figure(figObj);
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

% Figure setup
N_subfigs_per_fig = 24;
Nsubfigs = Nthetas * (Nthetas - 1) / 2;
Ncols = ceil( sqrt(2.5*min(Nsubfigs,N_subfigs_per_fig)) );
Nrows = ceil( min(Nsubfigs,N_subfigs_per_fig) / Ncols );
Nfigs = ceil(Nsubfigs / N_subfigs_per_fig);

% Create figure for parameter bivariate plots if it doesn't exist
for k = 1:Nfigs
    figObj = findobj('type','Figure','Name',['Bivariate Scatters ',num2str(k)]);
    if isempty(figObj)
        figs{k} = figure('units','Normalized','OuterPosition',[0 0 1 1], 'Name', ['Bivariate Scatters ',num2str(k)]);
    else
        figs{k} = figObj;
    end
end

% Initialise loop variables
c = 0;
cur_fig = 1;
% Move to current figure and clear it, then prepare axes
figure(figs{1});
clf;
ax = createAxes( Nsubfigs, Nrows, Ncols, [], true );

% Loop over parameter combinations
for i = 1:Nthetas-1
    for j = i+1:Nthetas

        % Increment counter
        c = c + 1;
        
        % Move to new figure if counter exceeds limit
        if c > N_subfigs_per_fig
            c = 1;
            cur_fig = cur_fig + 1;
            figure(figs{cur_fig});
            clf;
            ax = createAxes( Nsubfigs, Nrows, Ncols, [], true );
        end
        
        % Plot scatter
        plot(ax{c}, thetas(:,i), thetas(:,j), 'k.', 'MarkerSize', 15 );
        xlabel(ax{c}, theta_names{i}, 'FontSize', 16);
        ylabel(ax{c}, theta_names{j}, 'FontSize', 16);
        set(ax{c}, 'FontSize', 16);
    end
end


end