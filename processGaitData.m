% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6710137
% Based on the approach described in above paper
clear,close all
data = load("tom_trial_6.txt");
% parse into sensor locations
% 2-10 left ear; 11-19 right ear; 20-28 chest; 29:37 pocket
lear_data = data(:, 3:11);
rear_data = data(:, 12:20);
chest_data = data(:, 21:29);
pocket_data = data(:, 30:38);
% line up the axes of the sensors to allow for direct comparison
[lear_data, rear_data, chest_data, pocket_data] = alignEarData(lear_data, rear_data, chest_data, pocket_data);

using_right_ear = false;
if using_right_ear
    % perform SSA on the x-(Anterior-Posterior) axis
    [PC_AP, RC_AP] = ssa(rear_data(:, 1), false);
    % perform SSA on the y-(Superior-Inferior) axis
    [PC_SI, RC_SI] = ssa(rear_data(:, 2), false);
    % perform SSA on the x-(Medio-Lateral) axis
    [PC_ML, RC_ML] = ssa(rear_data(:, 3), false);
else
    % perform SSA on the x-(Anterior-Posterior) axis
    [PC_AP, RC_AP] = ssa(lear_data(:, 1), false);
    % perform SSA on the y-(Superior-Inferior) axis
    [PC_SI, RC_SI] = ssa(lear_data(:, 2), false);
    % perform SSA on the x-(Medio-Lateral) axis
    [PC_ML, RC_ML] = ssa(lear_data(:, 3), false);
%     [PC_AP, RC_AP] = ssa(chest_data(:, 1), false);
%     [PC_SI, RC_SI] = ssa(chest_data(:, 2), false);
%     [PC_ML, RC_ML] = ssa(chest_data(:, 3), false);
end

% save the dominant oscillations
dom_osc_AP = sum(RC_AP(:,2:2),2);
dom_osc_SI = sum(RC_SI(:,1:2),2);
dom_osc_ML = sum(RC_ML(:,1:2),2);
% save the overall filtered signals
filtered_AP = sum(RC_AP(:,2:end),2);
filtered_SI = sum(RC_SI(:,1:7),2);
filtered_ML = sum(RC_ML(:,1:15),2);
fig_filtered_signals = figure('Name','Filtered AP, SI and ML Signals');
plot(filtered_AP)
hold on
plot(filtered_SI)
plot(filtered_ML)
legend('AP', 'SI', 'ML')

% stack the ML and AP axes for inspection
fig_ML_SI = figure('Name', 'ML and SI axes');
subplot(2,1,1)
plot(dom_osc_ML);
hold on
plot(filtered_ML)
title({'ML Axis Signal with Dominant Oscillation (blue),', 'estimated heel strikes from AP (red dotted lines)', 'and revised estimate based on AP and SI (green dotted lines)'})
ylabel('Accel /g')
subplot(2,1,2)
[max_peak_amp, max_peak_loc] = findpeaks(dom_osc_AP);
[min_peak_amp, min_peak_loc] = findpeaks(-dom_osc_AP);
plot(dom_osc_AP);
hold on
title({'AP Axis Signal with Main Oscillation (blue)', 'and Estimated Heel Strikes (red dotted lines)'})
ylabel('Accel / g')
xlabel('Samples')
plot(filtered_AP)
plot(max_peak_loc, max_peak_amp, 'x');
plot(min_peak_loc, -min_peak_amp, 'o');
xline([min_peak_loc], '--r')
subplot(2,1,1)
xline([min_peak_loc(1:2:end)], '--r')
xline([min_peak_loc(2:2:end)], '--r')

% multiply AP and SI signals and select lowest vals
fig_AP_times_SI = figure('Name','AP times SI with AP dominant oscillation superimposed');
tau = 8; % width of interval around local minima
heel_contacts = zeros(1, length(min_peak_loc));
for i = 1:length(min_peak_loc)
    if (min_peak_loc(i)+tau < length(filtered_SI)) && (min_peak_loc(i)-tau > 0) 
        [~, temp_min_loc] = min(filtered_AP(min_peak_loc(i)-tau : min_peak_loc(i)+tau) .* filtered_SI(min_peak_loc(i)-tau : min_peak_loc(i)+tau));
        heel_contacts(i) = temp_min_loc + (min_peak_loc(i)-tau);
    end
end
plot(filtered_AP.* filtered_SI);
hold on
xline(heel_contacts, '--g')
plot(dom_osc_AP);
figure(fig_ML_SI)
xline(heel_contacts, '--g')

% now determine right and left labelling
% take 3 local minima from AP and compare these points on the ML
% if mean val between first 2 points > that between 2nd 2 points, then RHS
% ML_mean_1 = mean(dom_osc_ML(min_peak_loc(2):min_peak_loc(3)));
% ML_mean_2 = mean(dom_osc_ML(min_peak_loc(3):min_peak_loc(4)));
ML_mean_1 = mean(dom_osc_ML(heel_contacts(2):heel_contacts(3)));
ML_mean_2 = mean(dom_osc_ML(heel_contacts(3):heel_contacts(4)));
% let true mean top row is RHS and hence false means it is LHS
right_foot_first = (ML_mean_1 > ML_mean_2);
fig_heel_contacts_visualisation = figure('Name', 'Heel Contacts Visualisation');
plot(filtered_ML)
hold on
plot(dom_osc_ML)
xline(heel_contacts(2:4))
yline([ML_mean_1, ML_mean_2], '--')

% extract the gait cycles depending on left or right foot first
if right_foot_first
%     [main_gait_cycle, X_c, V] = extractGaitCycle(dom_osc_ML, min_peak_loc(2:2:end));
    [main_gait_cycle, X_c, V] = extractGaitCycle(filtered_ML, heel_contacts(2:2:end));
    X_c_LHCs = heel_contacts(3:2:end) - heel_contacts(2:2:end-1); % this ensures the LHC lines up for each cycle
else
%     [main_gait_cycle, X_c, V] = extractGaitCycle(dom_osc_ML, min_peak_loc(1:2:end));
    [main_gait_cycle, X_c, V] = extractGaitCycle(filtered_ML, heel_contacts(3:2:end));
    X_c_LHCs = heel_contacts(4:2:end) - heel_contacts(3:2:end);
end
fig_main_cycle = figure('Name', 'Main Gait Cycle Maxima and Minima i.e. Heel Contacts');
[main_cycle_max_peak_amp, main_cycle_max_peak_loc] = findpeaks(main_gait_cycle);
[main_cycle_min_peak_amp, main_cycle_min_peak_loc] = findpeaks(-main_gait_cycle);
plot(main_gait_cycle);
hold on
plot(main_cycle_max_peak_loc, main_cycle_max_peak_amp, 'x');
plot(main_cycle_min_peak_loc, -main_cycle_min_peak_amp, 'o');

% use lowest common subsequence to find the similarity
tau2 = 0.75; % similarity threshold - cycle discarded if below this
match_point = cell(1, size(X_c,1));
for j = 1:size(X_c,1)
    s1 = (main_gait_cycle - mean(main_gait_cycle)) / std(main_gait_cycle);
    s2 = (X_c(j,:) - mean(X_c(j,:))) / std(X_c(j,:));
    epsilon = 0.3*min(std(s1), std(s2)); % represents the degree of matching in the space, set to half input std
    delta = round(0.15 * min(length(main_gait_cycle), length(X_c(j,:)))); % temporal matching
    % calculate and check if similarity reaches threshold
    [similarity, matched_out] = MatLabLCSS(s1, s2, delta, epsilon, j);
    if similarity > tau2
        match_point{j} = matched_out;
    end
end

% loop through X_c (the gait cycles) to show LHC and hence find toe-offs
X_c_resampled = cell(1,size(X_c,1));
LTO_X_c = cell(1,3);
RTO_X_c = cell(1,3);
for i = 1:size(X_c,1)-1
    loc_LHC = X_c_LHCs(i+1);
    % find left toe-offs (LTO)
    [~, LTO_main_cycle] = max(main_gait_cycle(1:loc_LHC));
    if sum((match_point{i}(:,1) == LTO_main_cycle), 2) == 1
        LTO_X_c{i} = match_point{i}((match_point{i}(:,1) == LTO_main_cycle)',:);
    else
        nearestLTOIndex = findNearestIndex(LTO_main_cycle, match_point{i}(:,1));
        LTO_X_c{i} = match_point{i}(nearestLTOIndex,:);
    end
    figure('Name', sprintf('Gait Cycle %d Labelled', i));
    plot(X_c(i,:))
    hold on
    xline(loc_LHC, '-g')
    xline(LTO_X_c{i}(1,2), '--g')
    % find right toe-offs (RTO)
    [~, RTO_main_cycle] = min(main_gait_cycle(loc_LHC:end));
    if sum(match_point{i}(:,1) == (RTO_main_cycle + loc_LHC)) == 1
        RTO_X_c{i} = match_point{i}((match_point{i}(:,1) == (RTO_main_cycle + loc_LHC)),:);
    else
        nearestRTOIndex = findNearestIndex(RTO_main_cycle, match_point{i}(:,1));
        RTO_X_c{i} = match_point{i}(nearestRTOIndex + loc_LHC,:);
    end
    xline(RTO_X_c{i}(1,2), '--r')
end

% set the template r = a_1*exp{(t-b_1)^2/2*z12} + a_2*exp{(t-b_2)^2/2*z22}
% a_1, b_1, z12 = 10, 40, 600 ; a_2, b_2, z22 = 10, 120, 600 ; t = 1:150
% (half of the cycles, q = 310)
% I don't understand this part and it doesn't work for me
% t = 1:50;
% r = 10 * exp((t-40).^2/2*600) + 10 * exp((t - 120).^2 / 2*600);

% find maxima and minima of main gait cycle
% main_gait_peak_loc = max(main_cycle_max_peak_loc(main_cycle_max_peak_loc < 0.5*length(main_gait_cycle)))
% main_gait_trough_loc = min(main_cycle_min_peak_loc(main_cycle_min_peak_loc > 0.5*length(main_gait_cycle)))

% plot all the gait events on a single figure
fig_all_events = figure('Name', 'All Gait Events Pictured for Each Axis');
axis_data = {filtered_AP, filtered_SI, filtered_ML};
axis_titles = {'AP Signal', 'SI Signal', 'ML Signal'};
for j = 1:3
    subplot(3,1,j)
    plot(axis_data{j})
    hold on
    if right_foot_first%~xor(right_foot_first, using_right_ear)
        % store the RHCs and LHCs locations for the gait cycles
        overall_RHC_loc = heel_contacts(4:2:end);
        overall_LHC_loc = heel_contacts(3:2:end);
        % include first 2 on the plot
        xline(heel_contacts(1:2:end), '-g'); % left heel strikes
        xline(heel_contacts(2:2:end), '-r'); % right heel strikes
        % plot the LTOs
        for i = 1:size(X_c,1)-1
            LTO_locations(i) = LTO_X_c{i}(1,2) + heel_contacts(2*i);
            xline(LTO_X_c{i}(1,2) + heel_contacts(2*i), '--g')
        end
        % plot the RTOs
        for i = 1:size(X_c,1)-1
            RTO_locations(i) = RTO_X_c{i}(1,2) + heel_contacts(2*i);
            xline(RTO_X_c{i}(1,2) + heel_contacts(2*i), '--r')
        end
    else
        % store the RHCs and LHCs locations for the gait cycles
        overall_RHC_loc = heel_contacts(3:2:end);
        overall_LHC_loc = heel_contacts(4:2:end);
        % include first 2 here as well
        xline(heel_contacts(1:2:end), '-r'); % right heel strikes
        xline(heel_contacts(2:2:end), '-g'); % left heel strikes
         % plot the LTOs
        for i = 1:size(X_c,1)
            LTO_locations(i) = LTO_X_c{i}(1,2) + heel_contacts(2*i + 1);
            xline(LTO_X_c{i}(1,2) + heel_contacts(2*i + 1), '--g')
        end
        % plot the RTOs
        for i = 1:3 % size(X_c,1)
            RTO_locations(i) = RTO_X_c{i}(1,2) + heel_contacts(2*i + 1);
            xline(RTO_X_c{i}(1,2) + heel_contacts(2*i + 1), '--r')
        end
    end
    title(axis_titles{j});
    ylabel('Acceleration / g');
end
% figure(fig_all_events)
sgtitle('All Gait Events Pictured for Each Axis')

% calculate TSPs
if ~right_foot_first %xor(right_foot_first, using_right_ear)
    for cycle_num = 1:size(X_c,1)-1
        left_stance_time(cycle_num) = LTO_locations(cycle_num+1) - overall_LHC_loc(cycle_num);
        right_swing_time(cycle_num) = overall_RHC_loc(cycle_num+1) - RTO_locations(cycle_num);
    end
    for cycle_num = 1:size(X_c,1)
        left_swing_time(cycle_num) = overall_LHC_loc(cycle_num) - LTO_locations(cycle_num);
        right_stance_time(cycle_num) = RTO_locations(cycle_num) - overall_RHC_loc(cycle_num);
        left_stride_time(cycle_num) = overall_LHC_loc(cycle_num+1) - overall_LHC_loc(cycle_num);
        right_stride_time(cycle_num) = overall_RHC_loc(cycle_num+1) - overall_RHC_loc(cycle_num);
        step_asymmetry(cycle_num) = (overall_LHC_loc(cycle_num) - overall_RHC_loc(cycle_num)) / (overall_RHC_loc(cycle_num+1) - overall_LHC_loc(cycle_num));
    end
    % make dimensions agree
    left_stance_time = [0, left_stance_time];
    right_swing_time = [0, right_swing_time];
else
    for cycle_num = 1:size(X_c,1)-1
        left_stance_time(cycle_num) = LTO_locations(cycle_num+1) - overall_LHC_loc(cycle_num);
        left_swing_time(cycle_num) = overall_LHC_loc(cycle_num) - LTO_locations(cycle_num);
        right_stance_time(cycle_num) = RTO_locations(cycle_num+1) - overall_RHC_loc(cycle_num);
        right_stride_time(cycle_num) = overall_RHC_loc(cycle_num+1) - overall_RHC_loc(cycle_num);
    end
    for cycle_num = 1:size(X_c,1)
        right_swing_time(cycle_num) = overall_RHC_loc(cycle_num) - RTO_locations(cycle_num);
        left_stride_time(cycle_num) = overall_LHC_loc(cycle_num+1) - overall_LHC_loc(cycle_num);
        step_asymmetry(cycle_num) = (overall_RHC_loc(cycle_num) - overall_LHC_loc(cycle_num)) / (overall_LHC_loc(cycle_num+1) - overall_RHC_loc(cycle_num));
    end
    % make dimensions agree
    right_stance_time = [0, right_stance_time];
    right_stride_time = [0, right_stride_time];
    left_swing_time = [0, left_swing_time];
    left_stance_time = [0, left_stance_time];
end


% put them into a table
% TSPs = table(left_swing_time', right_swing_time', left_stance_time', right_stance_time', left_stride_time', right_stride_time',step_asymmetry',...
%      'VariableNames',{'LSwing','RSwing','LStance','RStance','LStride','RStride','Asymmetry'});

TSPs = table(left_stride_time', right_stride_time',step_asymmetry', 'VariableNames',{'LStride','RStride','Asymmetry'});

% save as a csv
writetable(TSPs, 'TSPs.csv', 'WriteRowNames', true);