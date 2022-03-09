function Biomarkers = computeBiomarkers(time,profile)

t=time;         %%%%% must be given in ms
V=profile;

% If the profile does not reach -30mV, we reject it as not an AP
% ALSO, if we have any complex numbers in our voltage, REJECT
AP_beginning = find(V>=-30,1);
if isempty(AP_beginning) || ~isreal(V) % Reject
    Class = 'RejectedNoDepolarisation';
    Biomarkers = nan(9,1);
else % Proceed to next check
    
    % Read out the peak locations and values
    [pks,locs]=findpeaks(V,'MinPeakProminence',0.1,'MinPeakHeight',-30);
    
    % Check now if there are more than two peaks - this is pretty much a
    % guarantee of a DAD, EAD or oscillations
    if length(pks) > 2
        Class = 'RejectedApparentDAD/EAD';
        Biomarkers = nan(9,1);
    else
        
        % Calculate the approximate derivative
        Vderiv = diff(V)./diff(t);
        
        % Find the flat part of the profile before the peak if possible
        [dVdtmin,min_loc] = min(abs(Vderiv(1:AP_beginning)));
        if dVdtmin < 0.005  % Sufficienctly flat location found
            RMP = V(min_loc);
        else
            RMP = min(V(1:AP_beginning));   % Sufficienctly flat location not found, RMP is minimal V pre-stimulus
        end
        
        % Now that the RMP has been calculated, the APA can also be calculated
        APA = pks(1) - RMP;
        
        
        %%% Now determine if the cell successfully repolarises
        repol_loc = find(V(AP_beginning:end) <= RMP+0.1*APA,1);
        
        if isempty(repol_loc)  % Failed to repolarise, reject
            Class = 'RejectedNoRepolarisation';
            Biomarkers = nan(9,1);
        else
            
            % Find the maximum upstroke velocity and its position
            [dVdtmax,max_loc] = max(Vderiv);
            t_init = t(max_loc);
            
            % Calculate the voltages at certain percentage repolarisations
            AP_20 = RMP + 0.8 * APA;
            AP_50 = RMP + 0.5 * APA;
            AP_90 = RMP + 0.1 * APA;
            
            % Find the locations at which the cell has repolarised to these extents
            repol_index_20 = find(V >= AP_20,1,'last');
            repol_index_50 = find(V >= AP_50,1,'last');
            repol_index_90 = find(V >= AP_90,1,'last');
            
            % Read out the time locations of these
            t_repol_20 = t(repol_index_20);
            t_repol_50 = t(repol_index_50);
            t_repol_90 = t(repol_index_90);
            
            % Now calculate APDs using the starting point of maximum upstroke
            % location
            APD_20 = t_repol_20 - t_init;
            APD_50 = t_repol_50 - t_init;
            APD_90 = t_repol_90 - t_init;
            
            % Find the time which is 20% of APD90
            t_20 = t_init + 0.2 * APD_90;
            t_20_loc = find(t >= t_20, 1);
            
            V20 = V(t_20_loc);
            
            % Now input biomarker values into the Biomarkers vector - the final two
            % regard dome peaks and so will only feature when a dome peak is
            % present
            if length(pks) == 1  % No dome peak (or dome peak less than -30mV which wasn't detected, who cares for now)
                Class = 'SinglePeak';
                Biomarkers = [APD_90; APD_50; APD_20; APA; RMP; dVdtmax; V20; NaN; NaN];
            else % Dome peak present, or a DAD
                
                % Check if this is a DAD or dome peak by testing for
                % repolarisation between peaks
                if repol_loc < locs(2) % Cell repolarises before second peak - DAD
                    Class = 'RejectedRepolarisesBeforeDomePeak';
                    Biomarkers = NaN * ones(9,1);
                else % Cell doesn't repolarise before second peak - NotchDome
                    
                    % Read out the voltage at the dome peak
                    V_domepeak = pks(2);
                    
                    % Find the minimum voltage value between the two peaks
                    V_interpeak_min = min(V(locs(1):locs(2)));
                    
                    Class = 'NotchDome';
                    Biomarkers = [APD_90; APD_50; APD_20; APA; RMP; dVdtmax; V20; V_interpeak_min; V_domepeak];
                end
            end
        end
    end
end

