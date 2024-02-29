
function [p_value] = MackSkill_SingleSession(Date)

%% Load the output structure

% Define where the excel spreadsheets are saved
Data_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\');

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
bsfr_morn_idxs = find(contains(xds_depth_excel.Properties.VariableNames, 'bsfr_morn'));
for ii = 1:length(bsfr_morn_idxs)
    % Find the unit
    unit_name = extractAfter(xds_depth_excel.Properties.VariableNames{bsfr_morn_idxs(ii)}, 'bsfr_morn_');
    % Morning baseline and movement phase firing rates
    bsfr_morn = mean(xds_depth_excel{:,bsfr_morn_idxs(ii)}, 'omitnan');
    mpfr_morn_idx = contains(xds_depth_excel.Properties.VariableNames, ...
        strcat('mpfr_morn_', unit_name));
    % Afternoon baseline and movement phase firing rates
    bsfr_noon_idx = contains(xds_depth_excel.Properties.VariableNames, ...
        strcat('bsfr_noon_', unit_name));
    bsfr_noon = mean(xds_depth_excel{:,bsfr_noon_idx}, 'omitnan');
    mpfr_noon_idx = contains(xds_depth_excel.Properties.VariableNames, ...
        strcat('mpfr_noon_', unit_name));
    % Calculate the depth of modulation
    depth_morn = xds_depth_excel{:,mpfr_morn_idx} - bsfr_morn;
    depth_morn(isnan(depth_morn)) = [];
    depth_noon = xds_depth_excel{:,mpfr_noon_idx} - bsfr_noon;
    depth_noon(isnan(depth_noon)) = [];
    % Trials
    temp_morn_trials = 1:length(depth_morn);
    temp_noon_trials = (1:length(depth_noon)) + length(depth_morn);
    % Final columns
    if isequal(ii, 1)
        depth_Column = cat(1,depth_morn,depth_noon);
        Time_Column = zeros(length(depth_Column),1);
        Time_Column(1:length(depth_morn),1) = 1;
        Time_Column(length(depth_morn):(length(depth_morn) + length(depth_noon)),1) = 2;
        Unit_Column = zeros(length(depth_Column),1) + ii;
        Trial_Column = cat(1, temp_morn_trials', temp_noon_trials');
    else
        temp_depth = cat(1,depth_morn,depth_noon);
        depth_Column = cat(1,depth_Column, temp_depth);
        temp_Time = zeros(length(temp_depth),1);
        temp_Time(1:length(depth_morn),1) = 1;
        temp_Time(length(depth_morn):(length(depth_morn) + length(depth_noon)),1) = 2;
        Time_Column = cat(1,Time_Column, temp_Time);
        temp_Unit = zeros(length(temp_depth),1) + ii;
        Unit_Column = cat(1,Unit_Column, temp_Unit);
        temp_Trial = cat(1, temp_morn_trials', temp_noon_trials');
        Trial_Column = cat(1, Trial_Column, temp_Trial);
    end
end

%% Stats

[p_value] = mackskill(depth_Column, Time_Column, Trial_Column);

p_value = p_value(1);



