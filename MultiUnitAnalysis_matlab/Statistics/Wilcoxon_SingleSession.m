
function [p_value] = Wilcoxon_SingleSession(Date)

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
depth_morn = zeros(length(depth_morn_idxs),1);
depth_noon = zeros(length(depth_morn_idxs),1);
for ii = 1:length(depth_morn_idxs)
    % Find the unit
    unit_name = extractAfter(xds_depth_excel.Properties.VariableNames{depth_morn_idxs(ii)}, ...
        'depth_morn_');
    % Morning depth of modulation
    depth_morn(ii,1) = mean(xds_depth_excel{:,depth_morn_idxs(ii)}, 'omitnan');
    % Afternoon depth of modulation
    depth_noon_idx = contains(xds_depth_excel.Properties.VariableNames, ...
        strcat('depth_noon_', unit_name));
    depth_noon(ii,1) =  mean(xds_depth_excel{:,depth_noon_idx}, 'omitnan');
end

%% Stats

[p_value, ~, ~] = ranksum(depth_morn, depth_noon);




