
function [p_value] = TwoAnova_SingleSession(Date)

%% Load the output structure

% Define where the excel spreadsheets are saved
Data_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\New\');

% Identify all the excel files in the data path
Excel_Path = strcat(Data_Path, '*.xlsx');
Excel_Files = dir(Excel_Path);
% Convert the struc to a table
Excel_Files_In_Path = struct2table(Excel_Files(~([Excel_Files.isdir])));

% Find the excel file
Excel_Choice = intersect(find(contains(Excel_Files_In_Path.name, Date)), ...
    find(~contains(Excel_Files_In_Path.name, '~')));

% Load the output table
table_path = strcat(Data_Path, Excel_Files_In_Path.name(Excel_Choice(1)));
xds_depth_excel = readtable(char(table_path));

%% Assign the variables
depth_morn_idxs = find(contains(xds_depth_excel.Properties.VariableNames, 'depth_morn'));
for ii = 1:length(depth_morn_idxs)
    % Find the unit
    unit_name = extractAfter(xds_depth_excel.Properties.VariableNames{depth_morn_idxs(ii)}, 'depth_morn_');
    % Morning depth of modulation
    depth_morn = xds_depth_excel{:,depth_morn_idxs(ii)};
    % Afternoon baseline and movement phase firing rates
    depth_noon_idx = contains(xds_depth_excel.Properties.VariableNames, ...
        strcat('depth_noon_', unit_name));
    depth_noon = xds_depth_excel{:,depth_noon_idx};
    % Calculate the depth of modulation
    depth_morn(isnan(depth_morn)) = [];
    depth_noon(isnan(depth_noon)) = [];
    % Trials
    temp_morn_trials = 1:length(depth_morn);
    temp_noon_trials = 1:length(depth_noon);
    % Final columns
    if isequal(ii, 1)
        depth_Column = cat(1,depth_morn,depth_noon);
        Time_Column = strings;
        Time_Column(1:length(depth_morn),1) = 'Morning';
        Time_Column(length(depth_morn)+1:length(depth_Column),1) = 'Afternoon';
        Trial_Column = cat(1, string(temp_morn_trials'), string(temp_noon_trials'));
    else
        temp_depth = cat(1,depth_morn,depth_noon);
        depth_Column = cat(1,depth_Column, temp_depth);
        temp_Time = strings;
        temp_Time(1:length(depth_morn),1) = 'Morning';
        temp_Time(length(depth_morn)+1:length(temp_depth),1) = 'Afternoon';
        Time_Column = cat(1,Time_Column, temp_Time);
        temp_Trial = cat(1, string(temp_morn_trials'), string(temp_noon_trials'));
        Trial_Column = cat(1, Trial_Column, temp_Trial);
    end
end

%% Stats
Trial_Column = char(Trial_Column);
Time_Column = char(Time_Column);
[rank_depth_Column, ~] = tiedrank(depth_Column);

p_value = anovan(rank_depth_Column, {Trial_Column Time_Column}, 'model', 2, ...
        'varnames', {'Trial', 'Time'});

p_value = p_value(1);


