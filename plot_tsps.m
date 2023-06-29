% plot the TSPs
clear, close all;

data = readtable("C:\Users\teri-\Downloads\tom-measurements\IMU TSPs\Stride Time.xlsx");

trials = data.Trial;
lstridemean = data.LStride_Mean;
rstridemean = data.RStride_Mean;
lstrideerrorl = data.LStride_err_l;
rstrideerrorl = data.RStride_err_l;
lstrideerrorh = data.LStride_err_h;
rstrideerrorh = data.RStride_err_h;

viconL = data.ViconL;
viconR = data.ViconR;



figure("Name", "Left Stride Times")
% scatter(trials(1:2:end), lstridemean(1:2:end), "green", "+")
hold on
rearplot = errorbar(trials(1:2:end), lstridemean(1:2:end), (lstridemean(1:2:end)-lstrideerrorl(1:2:end)), (lstrideerrorh(1:2:end)-lstridemean(1:2:end)),"xg");
learplot = errorbar(trials(2:2:end), lstridemean(2:2:end), (lstridemean(2:2:end)-lstrideerrorl(2:2:end)), (lstrideerrorh(2:2:end)-lstridemean(2:2:end)),"+r");
plot(trials(1:2:end), viconL(1:2:end), "ob")
xlabel("Trial Number")
ylabel("Step Time / centiseconds")
legend("Right Ear", "Left Ear", "Vicon","Location","north")
title({"Comparison of Calculated Left Step Times", " from Left and Right Ears with Vicon"})

figure("Name", "Right Stride Times")
hold on
rearplot = errorbar(trials(1:2:end), rstridemean(1:2:end), (rstridemean(1:2:end)-rstrideerrorl(1:2:end)), (rstrideerrorh(1:2:end)-rstridemean(1:2:end)),"xg");
learplot = errorbar(trials(2:2:end), rstridemean(2:2:end), (rstridemean(2:2:end)-rstrideerrorl(2:2:end)), (rstrideerrorh(2:2:end)-rstridemean(2:2:end)),"+r");
plot(trials(1:2:end), viconR(1:2:end), "ob")
xlabel("Trial Number")
ylabel("Step Time / centiseconds")
legend("Right Ear", "Left Ear", "Vicon","Location","north")
title({"Comparison of Calculated Right Step Times", " from Left and Right Ears with Vicon"})


swing_data = readtable("C:\Users\teri-\Downloads\tom-measurements\IMU TSPs\Swing Time.xlsx");

trials = swing_data.Trial;
lstridemean = swing_data.LSwing_Mean;
rstridemean = swing_data.RSwing_Mean;
lstrideerrorl = swing_data.LSwing_err_l;
rstrideerrorl = swing_data.RSwing_err_l;
lstrideerrorh = swing_data.LSwing_err_h;
rstrideerrorh = swing_data.RSwing_err_h;

viconL = swing_data.ViconR_LeftSingleSupport_;
viconR = swing_data.ViconL_RightSingleSupport_;



figure("Name", "Left Swing Times")
% scatter(trials(1:2:end), lstridemean(1:2:end), "green", "+")
hold on
rearplot = errorbar(trials(1:2:end), lstridemean(1:2:end), (lstridemean(1:2:end)-lstrideerrorl(1:2:end)), (lstrideerrorh(1:2:end)-lstridemean(1:2:end)),"xg");
learplot = errorbar(trials(2:2:end), lstridemean(2:2:end), (lstridemean(2:2:end)-lstrideerrorl(2:2:end)), (lstrideerrorh(2:2:end)-lstridemean(2:2:end)),"+r");
plot(trials(1:2:end), viconL(1:2:end), "ob")
xlabel("Trial Number")
ylabel("Swing Time / centiseconds")
legend("Right Ear", "Left Ear", "Vicon","Location","north")
title({"Comparison of Calculated Left Swing Times", " from Left and Right Ears with Vicon"})

figure("Name", "Right Swing Times")
hold on
rearplot = errorbar(trials(1:2:end), rstridemean(1:2:end), (rstridemean(1:2:end)-rstrideerrorl(1:2:end)), (rstrideerrorh(1:2:end)-rstridemean(1:2:end)),"xg");
learplot = errorbar(trials(2:2:end), rstridemean(2:2:end), (rstridemean(2:2:end)-rstrideerrorl(2:2:end)), (rstrideerrorh(2:2:end)-rstridemean(2:2:end)),"+r");
plot(trials(1:2:end), viconR(1:2:end), "ob")
xlabel("Trial Number")
ylabel("Swing Time / centiseconds")
legend("Right Ear", "Left Ear", "Vicon","Location","northwest")
title({"Comparison of Calculated Right Swing Times", " from Left and Right Ears with Vicon"})