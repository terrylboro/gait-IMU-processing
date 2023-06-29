% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6710137
% Based on the approach described in above paper
clear,close all
for trial = 5:5
    clearvars -except trial
    close all
%     data = load("tom_trial_5.txt");
    file_to_load = sprintf("tom_trial_%d.txt",trial);
    data = load(file_to_load);
    data = data(:, 3:11); % extract right ear data
    data(:, [3,6,9]) = -data(:, [3,6,9]); % lineup ML axis to fit convention
    data(:, [1,4,7]) = -data(:, [1,4,7]);
    
    % perform SSA on the x-(Anterior-Posterior),
    % y-(Superior-Inferior) and z-(Medio-Lateral) axes
    [PC_AP, RC_AP] = ssa(data(:, 1), false);
    [PC_SI, RC_SI] = ssa(data(:, 2), false);
    [PC_ML, RC_ML] = ssa(data(:, 3), false);
    
    % save the dominant oscillations
    dom_osc_AP = sum(RC_AP(:,2:2),2); 
    dom_osc_SI = sum(RC_SI(:,1:2),2);
    dom_osc_ML = sum(RC_ML(:,1:2),2);
    % save the overall filtered signals
    filtered_AP = sum(RC_AP(:,2:end),2); % osc 1 is trend we want rid of
    filtered_SI = sum(RC_SI(:,1:7),2); % these ranges are empirical
    filtered_ML = sum(RC_ML(:,1:15),2);
    % find the maxima and minima in the AP signal
    % to be used to identify foot contacts
    [AP_max_peak_amp, AP_max_peak_loc] = findpeaks(dom_osc_AP);
    [AP_min_peak_amp, AP_min_peak_loc] = findpeaks(-dom_osc_AP);
    
    % plot the processed signals
    fig_filtered_signals = figure('Name', 'Filtered AP, SI and ML Signals');
    % plot AP signal and dominant oscillation, with minima
    subplot(3,1,1); plot(filtered_AP); hold on; plot(dom_osc_AP);
    xline(AP_min_peak_loc, '--r'); % lines for foot contact est
    title({'AP Signal (Blue) with Dominant Oscillation (Solid Red)', 'and Estimated Foot Contacts (Dotted Red'}); ylabel('Accel /g');
    % plot SI filtered signal and dominant oscillation
    subplot(3,1,2); plot(filtered_SI); hold on; plot(dom_osc_SI);
    xline(AP_min_peak_loc, '--r'); % lines for foot contact est
    title({'SI Signal (Blue) with Dominant Oscillation (Solid Red)', 'and Estimated Foot Contacts (Dotted Red'}); ylabel('Accel /g');
    % plot the ML signals and the AP minima
    subplot(3,1,3); plot(filtered_ML); hold on;
    xline(AP_min_peak_loc, '--r'); % lines for foot contact est
    title('ML'); ylabel('Accel /g'); xlabel('Samples');
    sgtitle('Filtered AP, SI and ML Signals');
    
    % multiply AP and SI signals and select lowest vals
    fig_AP_times_SI = figure('Name','AP times SI with AP dominant oscillation superimposed');
    tau = 3; % width of interval around local minima
%     foot_contacts = zeros(1, length(AP_min_peak_loc));
    for i = 1:length(AP_min_peak_loc)
        if (AP_min_peak_loc(i)+tau < length(filtered_SI)) && (AP_min_peak_loc(i)-tau > 0) 
            [~, temp_min_loc] = min(filtered_AP(AP_min_peak_loc(i)-tau : AP_min_peak_loc(i)+tau) .* filtered_SI(AP_min_peak_loc(i)-tau : AP_min_peak_loc(i)+tau));
            foot_contacts(i) = temp_min_loc + (AP_min_peak_loc(i)-tau);
        end
    end
    plot(filtered_AP.* filtered_SI);
    hold on
    for i = 1:length(foot_contacts)
        str = string(foot_contacts(i));
        text(foot_contacts(i), -4.7, str)
    end
    if (foot_contacts(1) == 0), foot_contacts(1) = 2; 
    end
    
    % determine whether left or right foot first
    ML_mean_1 = mean(dom_osc_ML(foot_contacts(1):foot_contacts(2)));
    ML_mean_2 = mean(dom_osc_ML(foot_contacts(2):foot_contacts(3)));
    right_foot_first = (ML_mean_1 > ML_mean_2);
    % recolour the foot contacts on the graph (red = L, green = R)
    figure(fig_filtered_signals); subplot(3,1,3);
    if(right_foot_first)
        right_foot_contacts = foot_contacts(1:2:end);
        left_foot_contacts = foot_contacts(2:2:end);
        xline(right_foot_contacts, 'g');
    else
        right_foot_contacts = foot_contacts(2:2:end);
        left_foot_contacts = foot_contacts(1:2:end);
        xline(right_foot_contacts, 'g');
    end
    
    % extract gait cycles
    [main_gait_cycle, X_c, V] = extractGaitCycle(filtered_ML, right_foot_contacts);
    fig_main_cycle = figure('Name', 'Main Gait Cycle Maxima and Minima i.e. Heel Contacts');
    [main_cycle_max_peak_amp, main_cycle_max_peak_loc] = findpeaks(main_gait_cycle);
    [main_cycle_min_peak_amp, main_cycle_min_peak_loc] = findpeaks(-main_gait_cycle);
    plot(main_gait_cycle);
    hold on
    plot(main_cycle_max_peak_loc, main_cycle_max_peak_amp, 'x');
    plot(main_cycle_min_peak_loc, -main_cycle_min_peak_amp, 'o');
    
    % use lowest common subsequence to find the similarity
    tau2 = 0.75; % similarity threshold - cycle discarded if below this
    for j = 1:size(X_c,1)
        s1 = (main_gait_cycle - mean(main_gait_cycle)) / std(main_gait_cycle);
        s2 = (X_c(j,:) - mean(X_c(j,:))) / std(X_c(j,:));
        epsilon = 0.3*min(std(s1), std(s2)); % represents the degree of matching in the space, set to half input std
        delta = round(0.15 * min(length(main_gait_cycle), length(X_c(j,:)))); % temporal matching
        % calculate and check if similarity reaches threshold
        [similarity, matched_out] = MatLabLCSS(s1, s2, delta, epsilon, j, true);
        if similarity > tau2
            match_point{j} = matched_out;
        end
    end
    
    % loop through X_c (the gait cycles) to show LHC and hence find toe-offs
    X_c_resampled = cell(1,size(X_c,1));LTO_X_c = cell(1,3);RTO_X_c = cell(1,3);
    if right_foot_first
        % ensure L and R are aligned correctly then calculate the LHCs
        if length(right_foot_contacts) > length(left_foot_contacts)
            loc_LHC = left_foot_contacts - right_foot_contacts(1:end-1);
        else
            loc_LHC = left_foot_contacts - right_foot_contacts;
        end
    else
        if length(left_foot_contacts) > length(right_foot_contacts)
            loc_LHC = left_foot_contacts(2:end) - right_foot_contacts;
        else
            loc_LHC = left_foot_contacts(2:end) - right_foot_contacts(1:end-1);
        end
    end
    for i = 1:size(X_c,1)
        if (size(match_point{i}) > 1) % this ignores bad cycles
            % find left toe-offs (LTO)
            [~, LTO_main_cycle] = max(main_gait_cycle(1:loc_LHC(i)));
            if sum((match_point{i}(:,1) == LTO_main_cycle), 2) == 1
                LTO_X_c{i} = match_point{i}((match_point{i}(:,1) == LTO_main_cycle)',:);
            else % if there is no direct mapping, find closest point
                nearestLTOIndex = findNearestIndex(LTO_main_cycle, match_point{i}(:,1));
                LTO_X_c{i} = match_point{i}(nearestLTOIndex,:);
            end
            % find right toe-offs (RTO)
            [~, RTO_main_cycle] = min(main_gait_cycle(loc_LHC(i):end));
            if sum(match_point{i}(:,1) == (RTO_main_cycle + loc_LHC(i))) == 1
                RTO_X_c{i} = match_point{i}((match_point{i}(:,1) == (RTO_main_cycle + loc_LHC(i))),:);
            else % if there is no direct mapping, find closest point
                nearestRTOIndex = findNearestIndex(RTO_main_cycle, match_point{i}(:,1));
                RTO_X_c{i} = match_point{i}(nearestRTOIndex,:); %  + loc_LHC,
            end
            % plot the gait cycle and information
            figure('Name', sprintf('Gait Cycle %d Labelled', i));
            plot(X_c(i,:));hold on;xline(loc_LHC(i), '-r'); % left foot contact
            xline(LTO_X_c{i}(1,2), '--r');xline(RTO_X_c{i}(1,2), '--g') % toe-offs
        end
    end
    
    % plot all the gait events on a single figure
    fig_all_events = figure('Name', 'All Gait Events Pictured for Each Axis');
    axis_data = {filtered_AP, filtered_SI, filtered_ML};axis_titles = {'AP Signal', 'SI Signal', 'ML Signal'};
    for j = 1:3
        subplot(3,1,j)
        plot(axis_data{j})
        hold on
        xline(right_foot_contacts, '-g');
        xline(left_foot_contacts, '-r');
        % plot the LTOs
        for i = 1:size(X_c,1)
            if (size(match_point{i}) > 1) % this ignores bad cycles
                LTO_loc(i) = LTO_X_c{i}(1,2) + right_foot_contacts(i);
                xline(LTO_loc(i), '--r')
                RTO_loc(i) = RTO_X_c{i}(1,2) + right_foot_contacts(i);
                xline(RTO_loc(i), '--g')
            end
        end
        title(axis_titles{j}); ylabel('Acceleration / g');
    end
    sgtitle('All Gait Events Pictured for Each Axis')
    
    TSPs = calculateTSPs(left_foot_contacts, right_foot_contacts, LTO_loc, RTO_loc, right_foot_first);
    tsps_filename = sprintf("lear-tsps-%d.xlsx",trial);
    writetable(TSPs, tsps_filename);
end

