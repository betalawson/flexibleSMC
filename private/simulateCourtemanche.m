function [t,y] = simulateCourtemanche(params)
% This function simulates the Courtemanche et al. atrial model using the
% parameters given in pStar, at the times given in t


%%% Pulsing parameters
settle_time = 2500;   % Amount of time to wait before pacing (in ms)
pulse_BCL = 1000;     % Amount of time between stimulus pulses
num_prepulses = 4;    % Number of pulses used before collecting output

%%% Stimulus parameters
stim_amp = 2210;
stim_dur = 2;

% Initialise the CRN model
y0 = courtemancheInitialise();

% First, run the model for the settle time (no stimulus)
[~,~,y_end] = simulateCourtemanchePulse(params, y0, 0:1:settle_time, 0, 0);

% Now run the pre-pulses
for k = 1:num_prepulses
    [~,~,y_end] = simulateCourtemanchePulse(params, y_end, 0:1:pulse_BCL, stim_amp, stim_dur);
end

% Now run the final stimulus and collect the output
[t,y,~] = simulateCourtemanchePulse(params, y_end, 0:1:500, stim_amp, stim_dur);

