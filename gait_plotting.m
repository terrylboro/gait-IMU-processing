clear, close all;
% load the trial data
data = load("tom_trial_4.txt");
% perform filtering
% data = lowpass(data,30,100);

% parse into sensor locations
% 2-10 left ear; 11-19 right ear; 20-28 chest; 29:37 pocket
lear_data = data(:, 3:11);
rear_data = data(:, 12:20);
chest_data = data(:, 21:29);
pocket_data = data(:, 30:38);
% line up the axes of the sensors to allow for direct comparison
[lear_data, rear_data, chest_data, pocket_data] = alignEarData(lear_data, rear_data, chest_data, pocket_data);

% create an array with all in so we can loop through it in future
data_array = {lear_data, rear_data, chest_data, pocket_data};

% Compute CWT using Morlet wavelet and scales from 1 to 100
% scales = 1:100;
% cwt(chest_data(:,2), scales ,'gaus1');

% cwt = cwtft(chest_data(:,2).', 'wavelet', 'dog', 'scales', scales);
% Compute differentiation of CWT coefficients
% dcwt = diff(abs(cwt.cfs), 1, 2);
% % Find maxima and minima of differentiated coefficients
% maxima = findpeaks(dcwt);
% minima = findpeaks(-dcwt);

% plot accelerometer vals from all sensors
plotIMUCols('X-Accel', data_array, 1);
plotIMUCols('Y-Accel', data_array, 2);
plotIMUCols('Z-Accel', data_array, 3);

% plot gyro vals from all sensors
plotIMUCols('X-Gyr', data_array, 4);
plotIMUCols('Y-Gyr', data_array, 5);
plotIMUCols('Z-Gyr', data_array, 6);

% plot mag vals from all sensors
plotIMUCols('X-Mag', data_array, 7);
plotIMUCols('Y-Mag', data_array, 8);
plotIMUCols('Z-Mag', data_array, 9);




