function ax = CreateAxes(N_figs, N_rows, N_cols, ax_properties, single_index)
% This function creates a new set of axes for plotting - axes created are
% equal in size and positioned according to the sizing properties provided.
% Properties that aren't provided are assigned default values, which may be
% changed here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default.margin = 0.05;
default.xgap = 0.05;
default.ygap = 0.1;
default.leftspace = 0;
default.topspace = 0;
default.rightspace = 0;
default.bottomspace = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise axis properties to a dummy value if not provided
if nargin < 4
    ax_properties = [];
end
if nargin < 5
    single_index = false;
end

% Read out the requested properties or assign default values where they are
% not provided
if isfield(ax_properties,'margin')
    margin = ax_properties.margin;
else
    margin = default.margin;
end
if isfield(ax_properties,'xgap')
    xgap = ax_properties.xgap;
else
    xgap = default.xgap;
end
if isfield(ax_properties,'ygap')
    ygap = ax_properties.ygap;
else
    ygap = default.ygap;
end
if isfield(ax_properties,'leftspace')
    leftspace = ax_properties.leftspace;
else
    leftspace = default.leftspace;
end
if isfield(ax_properties,'topspace')
    topspace = ax_properties.topspace;
else
    topspace = default.topspace;
end
if isfield(ax_properties,'rightspace')
    rightspace = ax_properties.rightspace;
else
    rightspace = default.rightspace;
end
if isfield(ax_properties,'bottomspace')
    bottomspace = ax_properties.bottomspace;
else
    bottomspace = default.bottomspace;
end

% Calculate the width and height of each axes so that all will fit
dx = (1 - 2*margin - (N_cols-1)*xgap - leftspace - rightspace) / N_cols;
dy = (1 - 2*margin - (N_rows-1)*ygap - topspace - bottomspace) / N_rows;

% Create the axes objects
ax = cell(N_rows, N_cols);
k = 0;
for i = 1:N_rows
    for j = 1:N_cols
        
        % Find starting (x,y) co-ordinate for this axis
        xloc = margin + leftspace + (j-1)*(xgap+dx);           % Work from left and add
        yloc = 1 - margin - topspace - (i-1)*(ygap+dy) - dy;   % Work from top and subtract
        
        % Create the actual axes if this is within the requested number of
        % figures
        k = k+1;
        if k <= N_figs
            
            if single_index
                ax{k} = axes('Position',[xloc yloc dx dy]);
            else
                ax{i,j} = axes('Position',[xloc yloc dx dy]);
            end
        end

    end
end