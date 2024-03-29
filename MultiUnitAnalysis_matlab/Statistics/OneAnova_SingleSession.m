
function [p_value] = OneAnova_SingleSession(Date)

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
    % Afternoon depth of modulation
    depth_noon_idx = contains(xds_depth_excel.Properties.VariableNames, ...
        strcat('depth_noon_', unit_name));
    depth_noon = xds_depth_excel{:,depth_noon_idx};
    % Calculate the depth of modulation
    depth_morn(isnan(depth_morn)) = [];
    depth_noon(isnan(depth_noon)) = [];

    % Final columns
    if isequal(ii, 1)
        depth_Matrix = cat(1,depth_morn,depth_noon);
        Time_Column = strings;
        Time_Column(1:length(depth_morn),1) = 'Morning';
        Time_Column(length(depth_morn)+1:length(depth_Matrix),1) = 'Afternoon';
    else
        temp_depth = cat(1,depth_morn,depth_noon);
        depth_Matrix = cat(2,depth_Matrix, temp_depth);
        %temp_time = strings;
        %temp_time(1:length(depth_morn),1) = 'Morning';
        %temp_time(length(depth_morn)+1:length(temp_depth),1) = 'Afternoon';
        %Time_Column = cat(1,Time_Column, temp_time);
    end
end

%% Stats
[rank_depth_Matrix, ~] = tiedrank(depth_Matrix);

[p_value, ~, ~] = anova1(rank_depth_Matrix', Time_Column, 'off');

%figure
%hold on
%[~,~,~,~] = multcompare(stats, "Dimension",[1 2]);

