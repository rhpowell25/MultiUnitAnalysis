function Scatter_Compare_PerTrial(Date)

%% Load the output structure

% Define where the excel spreadsheets are saved
Data_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\');

% Identify all the excel files in the data path
Excel_Path = strcat(Data_Path, 'Old\*.xlsx');
Excel_Files = dir(Excel_Path);
% Convert the struc to a table
Excel_Files_In_Path = struct2table(Excel_Files(~([Excel_Files.isdir])));

% Find the excel file
Excel_Choice = intersect(find(contains(Excel_Files_In_Path.name, Date)), ...
    find(~contains(Excel_Files_In_Path.name, '~')));

% Load the output table
table_path = strcat(Data_Path, 'Old\', Excel_Files_In_Path.name(Excel_Choice(1)));
Old_depth_excel = readtable(char(table_path));

% Identify all the excel files in the data path
Excel_Path = strcat(Data_Path, 'New\*.xlsx');
Excel_Files = dir(Excel_Path);
% Convert the struc to a table
Excel_Files_In_Path = struct2table(Excel_Files(~([Excel_Files.isdir])));

% Find the excel file
Excel_Choice = intersect(find(contains(Excel_Files_In_Path.name, Date)), ...
    find(~contains(Excel_Files_In_Path.name, '~')));

% Load the output table
table_path = strcat(Data_Path, 'New\', Excel_Files_In_Path.name(Excel_Choice(1)));
New_depth_excel = readtable(char(table_path));

%% Assign the variables
depth_morn_idxs = find(contains(Old_depth_excel.Properties.VariableNames, 'depth_morn'));
for ii = 1:length(depth_morn_idxs)
    % Find the unit
    unit_name = extractAfter(Old_depth_excel.Properties.VariableNames{depth_morn_idxs(ii)}, ...
        'depth_morn_');
    % Morning depth of modulation
    depth_morn = Old_depth_excel{:,depth_morn_idxs(ii)};
    % Afternoon depth of modulation
    depth_noon_idx = contains(Old_depth_excel.Properties.VariableNames, ...
        strcat('depth_noon_', unit_name));
    depth_noon = Old_depth_excel{:,depth_noon_idx};
    % Remove any NaN's
    depth_morn(isnan(depth_morn)) = [];
    depth_noon(isnan(depth_noon)) = [];
    % Final columns
    if isequal(ii, 1)
        Old_depth_Column = cat(1,depth_morn,depth_noon);
        Old_Time_Column = zeros(length(Old_depth_Column),1);
        Old_Time_Column(1:length(depth_morn),1) = 1;
        Old_Time_Column(length(depth_morn):(length(depth_morn) + length(depth_noon)),1) = 2;
        Old_Unit_Column = strings;
        Old_Unit_Column(1:length(Old_depth_Column),1) = unit_name;
    else
        temp_depth = cat(1,depth_morn,depth_noon);
        Old_depth_Column = cat(1,Old_depth_Column, temp_depth);
        temp_Time = zeros(length(temp_depth),1);
        temp_Time(1:length(depth_morn),1) = 1;
        temp_Time(length(depth_morn):(length(depth_morn) + length(depth_noon)),1) = 2;
        Old_Time_Column = cat(1,Old_Time_Column, temp_Time);
        temp_Unit = strings;
        temp_Unit(1:length(temp_depth),1) = unit_name;
        Old_Unit_Column = cat(1,Old_Unit_Column, temp_Unit);
    end
end

depth_morn_idxs = find(contains(New_depth_excel.Properties.VariableNames, 'depth_morn'));
for ii = 1:length(depth_morn_idxs)
    % Find the unit
    unit_name = extractAfter(New_depth_excel.Properties.VariableNames{depth_morn_idxs(ii)}, ...
        'depth_morn_');
    % Morning depth of modulation
    depth_morn = New_depth_excel{:,depth_morn_idxs(ii)};
    % Afternoon depth of modulation
    depth_noon_idx = contains(New_depth_excel.Properties.VariableNames, ...
        strcat('depth_noon_', unit_name));
    depth_noon = New_depth_excel{:,depth_noon_idx};
    % Remove any NaN's
    depth_morn(isnan(depth_morn)) = [];
    depth_noon(isnan(depth_noon)) = [];
    % Final columns
    if isequal(ii, 1)
        New_depth_Column = cat(1,depth_morn, depth_noon);
        New_Time_Column = zeros(length(New_depth_Column),1);
        New_Time_Column(1:length(depth_morn),1) = 1;
        New_Time_Column(length(depth_morn):(length(depth_morn) + length(depth_noon)),1) = 2;
        New_Unit_Column = strings;
        New_Unit_Column(1:length(New_depth_Column),1) = unit_name;
    else
        temp_depth = cat(1,depth_morn,depth_noon);
        New_depth_Column = cat(1,New_depth_Column, temp_depth);
        temp_Time = zeros(length(temp_depth),1);
        temp_Time(1:length(depth_morn),1) = 1;
        temp_Time(length(depth_morn):(length(depth_morn) + length(depth_noon)),1) = 2;
        New_Time_Column = cat(1,New_Time_Column, temp_Time);
        temp_Unit = strings;
        temp_Unit(1:length(temp_depth),1) = unit_name;
        New_Unit_Column = cat(1,New_Unit_Column, temp_Unit);
    end
end

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? 
% ('Baseline', 'Peak', 'Depth', 'Depth_Change', 'alignment', ')?
fire_rate_phase = 'Depth';

% Do you want to look at the morning or afternoon ('Morn', 'Noon')
Morn_vs_Noon = 'Noon';

%% Reassign variables according to what you're plotting

if strcmp(fire_rate_phase, 'Depth')
    disp('Depth of Modulation')
    if strcmp(Morn_vs_Noon, 'Morn')
        fire_rate_first = Old_depth_Column(Old_Time_Column == 1);
        fire_rate_second = New_depth_Column(New_Time_Column == 1);
        Old_Unit_Column = Old_Unit_Column(Old_Time_Column == 1);
        New_Unit_Column = New_Unit_Column(New_Time_Column == 1);
    elseif strcmp(Morn_vs_Noon, 'Noon')
        fire_rate_first = Old_depth_Column(Old_Time_Column == 2);
        fire_rate_second = New_depth_Column(New_Time_Column == 2);
        Old_Unit_Column = Old_Unit_Column(Old_Time_Column == 2);
        New_Unit_Column = New_Unit_Column(New_Time_Column == 2);
    end
end

%% Some variable extraction & definitions

% Font & plotting specifications
[Plot_Params] = Plot_Specs;

axis_expansion = 0.01;

% Scatter Marker Shapes
marker_metric = '.';

% Scatter Marker Colors
color_metric = 0;

% Scatter Marker sizes
sz = 500;

%% Plot the Depth of Modulation Scatter

scatter_fig = figure;
scatter_fig.Position = [200 50 Plot_Params.fig_size Plot_Params.fig_size];
hold on

% Label the axis
if strcmp(fire_rate_phase, 'Depth')
    xlabel('First Depth of Modulation (Hz)', 'FontSize', Plot_Params.label_font_size);
    ylabel('Second Depth of Modulation (Hz)', 'FontSize', Plot_Params.label_font_size);
end

% Calculate the axis limits
round_digit = 1;
min_first = min(fire_rate_first);
min_second = min(fire_rate_second);
axis_min = round(min(min_first, min_second)/round_digit)*round_digit;
max_first = max(fire_rate_first);
max_second = max(fire_rate_second);
axis_max = round(max(max_first, max_second)/round_digit)*round_digit;

% Draw the unity line 
line([axis_min - axis_expansion, axis_max + axis_expansion],...
    [axis_min - axis_expansion, axis_max + axis_expansion], ... 
    'Color', 'k', 'Linewidth', 1, 'Linestyle','--')

conserved_units = intersect(Old_Unit_Column, New_Unit_Column);
for ii = 1:length(conserved_units)
    Old_Unit_idx = strcmp(Old_Unit_Column, conserved_units{ii});
    New_Unit_idx = strcmp(New_Unit_Column, conserved_units{ii});
    scatter(fire_rate_first(Old_Unit_idx), fire_rate_second(New_Unit_idx), ...
        sz, marker_metric, 'MarkerEdgeColor', [color_metric 0 0], ...
        'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);
end

% Display the means
fire_rate_mean_first = mean(fire_rate_first);
fire_rate_mean_second = mean(fire_rate_second);
fprintf('First average is %0.2f \n', fire_rate_mean_first)
fprintf('Second average is %0.2f \n', fire_rate_mean_second)

% Set the axis
xlim([axis_min - axis_expansion, axis_max + axis_expansion])
ylim([axis_min - axis_expansion, axis_max + axis_expansion])




