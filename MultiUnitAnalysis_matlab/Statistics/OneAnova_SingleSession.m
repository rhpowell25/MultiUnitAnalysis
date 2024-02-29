
function [p_value] = Anova_OneWay_SingleSession(Date)

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

    % Final columns
    if isequal(ii, 1)
        depth_Matrix = cat(1,depth_morn,depth_noon);
        Time_Column = strings;
        Time_Column(1:length(depth_morn),1) = 'Morning';
        Time_Column(length(depth_morn):(length(depth_morn) + length(depth_noon)),1) = 'Afternoon';
    else
        temp_depth = cat(1,depth_morn,depth_noon);
        depth_Matrix = cat(2,depth_Matrix, temp_depth);
    end
end

%% Stats
[rank_depth_Matrix, ~] = tiedrank(depth_Matrix);

depth_Matrix(strcmp(Time_Column, 'Morning'),:)

%mean
ttest2(depth_Matrix(strcmp(Time_Column, 'Morning'),:), depth_Matrix(strcmp(Time_Column, 'Afternoon'),:))
[~, tbl, stats] = anova1(rank_depth_Matrix', Time_Column);
%[p_value, ~, stats] = anova1(depth_Matrix', Time_Column);
%[~, tbl, stats] = kruskalwallis(depth_Matrix', Time_Column);
p_value = cell2mat(tbl(2,5));

figure
hold on
[~,~,~,~] = multcompare(stats, "Dimension",[1 2]);

