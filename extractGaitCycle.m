function [gait_cycle, X_c, V] = extractGaitCycle(ML_data, RHS_locs)
    % separate ML axis after previous trend removal and mean correction (9)
    X_c = [];%zeros(length(RHS_locs)-1);

       % q = max{RHC(i + j +1)−RHC(i + j)}j=0,...,N −1 all i+ 1
    diff_len = max(diff(RHS_locs));
    for loc = 1:length(RHS_locs)-1
        if (RHS_locs(loc) + diff_len) < length(ML_data)
            X_c = [X_c; ML_data(RHS_locs(loc):RHS_locs(loc)+diff_len)'];
        figure
        plot(ML_data(RHS_locs(loc):RHS_locs(loc)+diff_len)')
        end
    end
    % apply SVD (10)
    [~,~,V] = svd(X_c, "econ");
    % gc = v1 /norm(v1 ) (11)
    % where v1 is the first vector in V i.e. largest eigenvalue
    gait_cycle = V(:,1) / norm(V(:,1));
    % correct error where Matlab reports sign-inverted V
    % https://uk.mathworks.com/help/matlab/ref/double.svd.html;jsessionid=52ee6269c7db12fabcb0774a4306
    if sign(gait_cycle(7) - gait_cycle(3)) < 0
        gait_cycle = -gait_cycle;
    end
end