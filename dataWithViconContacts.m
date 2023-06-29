clear,close all
trials = [2,5,6,7,8];

contacts_matrix = [];
selectedEar = "Right";

for i = 1:length(trials)
    trial_num = trials(i);
    contacts_filename = "15_CU_A096391_01_0010\15_CU_A096391_01_00" + sprintf("%02d",trial_num) + ".csv";
    data_filename = "tom_trial_" + trial_num + ".txt";
    footContacts = readmatrix(contacts_filename, 'Range', 'D34:D50');
    contactSides = readcell(contacts_filename, 'Range', 'B34:B50');
    strikeOrOff = readcell(contacts_filename, 'Range', 'C34:C50');

    first_NaN_index = find(isnan(footContacts), 1);
    footContacts = footContacts(1:first_NaN_index-1);
    contactSides = contactSides(1:first_NaN_index-1);
    strikeOrOff = strikeOrOff(1:first_NaN_index-1);
    
    data = load(data_filename);

    r_data = data(:, 12:14); % extract right ear data
    r_data(:,2) = -r_data(:,2)
    r_data(:, 3) = -r_data(:, 3); % lineup ML axis to fit convention
    l_data = data(:, 3:5); % extract left ear data
    l_data(:,2) = -l_data(:,2)
    l_data(:, 1) = -l_data(:, 1); % lineup AP axis to fit convention

    rfig = figure('Name', 'All Gait Events Pictured for Each Axis - Right');
    r_axis_data = {r_data(:,1),r_data(:,2),r_data(:,3)};l_axis_data = {l_data(:,1),l_data(:,2),l_data(:,3)};
    r_axis_titles = {'Right AP Signal', 'Right SI Signal', 'Right ML Signal'};
    l_axis_titles = {'Left AP Signal', 'Left SI Signal', 'Left ML Signal'};
    for j = 1:3
        subplot(3,1,j)
        plot(r_axis_data{j})
        hold on
        xline(footContacts(matches(contactSides,"Left") & matches(strikeOrOff,"Foot Strike")) * 100, 'r');
        xline(footContacts(matches(contactSides,"Left") & matches(strikeOrOff,"Foot Off")) * 100, '--r');
        xline(footContacts(matches(contactSides,"Right") & matches(strikeOrOff,"Foot Strike")) * 100, 'g');
        xline(footContacts(matches(contactSides,"Right") & matches(strikeOrOff,"Foot Off")) * 100, '--g');
        title(r_axis_titles{j}); ylabel('Acceleration / g');
    end
    xlabel('Samples')
    sgtitle({'All Gait Events Pictured for Each Axis', "Right Ear Trial "+string(trial_num)})
    saveas(rfig,"Right Ear Trial "+string(trial_num)+".png")

    lfig = figure('Name', 'All Gait Events Pictured for Each Axis - Left');
    for k = 1:3
        subplot(3,1,k)
        plot(l_axis_data{k})
        hold on
        xline(footContacts(matches(contactSides,"Left") & matches(strikeOrOff,"Foot Strike")) * 100, 'r');
        xline(footContacts(matches(contactSides,"Left") & matches(strikeOrOff,"Foot Off")) * 100, '--r');
        xline(footContacts(matches(contactSides,"Right") & matches(strikeOrOff,"Foot Strike")) * 100, 'g');
        xline(footContacts(matches(contactSides,"Right") & matches(strikeOrOff,"Foot Off")) * 100, '--g');
        title(l_axis_titles{k}); ylabel('Acceleration / g');
    end
    xlabel('Samples')
    sgtitle({'All Gait Events Pictured for Each Axis,', "Left Ear Trial "+string(trial_num)})
    saveas(lfig,"Left Ear Trial "+string(trial_num)+".png")

    contacts_matrix = [contacts_matrix; footContacts];
end