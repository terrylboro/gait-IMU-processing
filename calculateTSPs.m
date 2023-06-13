function [TSPs] = calculateTSPs(LHC, RHC, LTO, RTO, right_foot_first)
    L_stride_time = diff(LHC);
    R_stride_time = diff(RHC);
    % calculations will depend on number of gait cycles
    % which is reflected in length of LTO/RTO
    if right_foot_first
        L_stance_time = LTO(2:end) - LHC(LHC < max(LTO));
        R_stance_time = RTO - RHC(RHC < max(RTO));
        L_swing_time = LHC(1:length(LTO)) - LTO;
        R_swing_time = RHC(2:1+length(RTO)) - RTO;
    %     step_asymmetry = R_stride_time / (overall_LHC_loc(cycle_num+1) - overall_RHC_loc(cycle_num));
        % add padding
    else
        L_stance_time = LTO - LHC(1:length(LTO));
        R_stance_time = RTO - RHC(RHC < max(RTO));
        L_swing_time = LHC(2:1+length(LTO)) - LTO;
        R_swing_time = RHC(2:1+length(RTO)) - RTO;
    end
    
    % perform padding to allow publishing in a table
    TSPs_cells = {L_swing_time, R_swing_time, L_stance_time, R_stance_time, L_stride_time, R_stride_time};
    TSP_variable_lengths = cellfun('size',TSPs_cells,2);
    max_TSP_variable_length = max(TSP_variable_lengths);
    for i = 1:length(TSPs_cells)
        if TSP_variable_lengths(i) < max_TSP_variable_length
            TSPs_cells{i} = [0, TSPs_cells{i}];
        end
    end
    TSPs = cell2table(TSPs_cells, "VariableNames",{'LSwing','RSwing','LStance','RStance','LStride','RStride'});
end