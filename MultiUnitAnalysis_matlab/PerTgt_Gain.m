%% Load the output structures
clear
clc

% What Drug Do You Want To Load? ('Caff', 'Cyp', 'Lex', 'Con') or ('YYYYMMDD')
Drug_Choice = 'All';
% What targets do you want to use ('Max' vs 'Min')
Tgt_Choice = 'All';

[~, ~, depth, ~, ~, ~, ~, ~, ~, ~, ~, ~, pref_dir, tgt_center, ~, ...
    file_name, all_unit_names] = Fast_Load_Depth(Drug_Choice, Tgt_Choice);

%% Some of the plotting specifications
% Individual Sessions or All Sessions? ('per_trial' vs 'all_trials')
trial_sessions = 'per_trial';

% Which task do you want to plot ('PG' vs. 'WS', vs 'Both)?
trial_task = 'PG';

% What preferred direction do you want to plot(-90, 0, 90, 180, 'All')
plot_dir = 'All';

% What minimum depth of modulation do you want to observe (#, or NaN)
depth_min = 0;

% Do you want to plot a line of best fit? (1 = Yes or 0 = No)
error_choice = 1;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
Save_Figs = 0;
if ~isequal(Save_Figs, 0)
    close all
end

%% Some variable extraction & definitions

% Font specifications
label_font_size = 15;
title_font_size = 12;
font_name = 'Arial';
fig_size_x = 400;
fig_size_y = 400;

% Save Counter
ss = 1;

% Scatter Marker sizes
marker_size = 250;
% Scatter Marker Shapes
marker_shape = '.';

%% Only use trial tasks selected
if strcmp(trial_task, 'PG')
    disp('Powergrasp Only')
    powergrasp_idx = find(~contains(file_name, 'PG'));
    all_unit_names(powergrasp_idx) = [];
    depth(powergrasp_idx) = [];
    pref_dir(powergrasp_idx) = [];
    tgt_center(powergrasp_idx) = [];
    file_name(powergrasp_idx) = [];
end

if strcmp(trial_task, 'WS')
    disp('Wrist w/ Spring Only')
    wrist_idx = find(~contains(file_name, 'WS'));
    all_unit_names(wrist_idx) = [];
    depth(wrist_idx) = [];
    pref_dir(wrist_idx) = [];
    tgt_center(wrist_idx) = [];
    file_name(wrist_idx) = [];
end

%% Only use the preferred direction selected

if ~strcmp(plot_dir, 'All')
    fprintf('Preferred Direciton of %.0f° Only \n', plot_dir)
    for ii = 1:length(pref_dir)
        pref_dir_idx = find(pref_dir{ii,1} ~= plot_dir);
        all_unit_names{ii,1}(pref_dir_idx) = [];
        depth{ii,1}(pref_dir_idx) = [];
        tgt_center{ii,1}(pref_dir_idx) = [];
    end
end

%% Remove any units that fall below the minimum selected depth of modulation

if ~isnan(depth_min)
    min_depth_violations = struct([]);
    for ii = 1:length(all_unit_names)
        min_depth_violations{ii,1} = find(depth{ii,1} < depth_min);
        all_unit_names{ii,1}(min_depth_violations{ii,1}) = [];
        depth{ii,1}(min_depth_violations{ii,1}) = [];
        tgt_center{ii,1}(min_depth_violations{ii,1}) = [];
    end
end

%% Merge all the information (Only if you selected 'all_trial')
if strcmp(trial_sessions, 'all_trials')
    % Find the amount of units each experiment uses
    units_per_experiment = zeros(length(tgt_center),1);
    for ii = 1:length(units_per_experiment)
        units_per_experiment(ii) = length(tgt_center{ii});
    end
    total_units = sum(units_per_experiment);

    % Concatenate all the information
    all_trials_depth = zeros(total_units, 1);
    all_trials_targets = zeros(total_units, 1);
    cc = 1;
    for xx = 1:length(tgt_center)
        for jj = 1:length(tgt_center{xx})
            all_trials_depth(cc,1) = depth{xx}(jj);
            all_trials_targets(cc,1) = tgt_center{xx}(jj);
            cc = cc + 1;
        end
    end
    % Rename the now merged variables
    depth = struct([]);
    depth{1,1} = all_trials_depth;
    tgt_center = struct([]);
    tgt_center{1,1} = all_trials_targets;
end

%% Find the mean modulation at each target center
avg_depth = struct([]);
target_centers = struct([]);
for ii = 1:length(depth)
    target_centers{ii,1} = unique(tgt_center{ii,1});
    for tt = 1:length(target_centers{ii,1})
        tgt_idx = find(tgt_center{ii,1} == target_centers{ii,1}(tt));
        avg_depth{ii,1}(tt,1) = mean(depth{ii,1}(tgt_idx));
    end
end

%% Plot the depth of modulation as a function of target center

scatter_fig = figure;
scatter_fig.Position = [200 50 fig_size_x fig_size_y];
hold on

% Set the title
title_string = strcat('Mean Modulation', {', '}, trial_task);
title(title_string, 'FontSize', title_font_size)
% Label the axis
xlabel('Target Distance', 'FontSize', label_font_size + 5);
ylabel('Depth of Modulation (Hz)', 'FontSize', label_font_size + 5);

% Axis Editing
figure_axes = gca;
% Set ticks to outside
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off')
% Set the tick label font size
figure_axes.FontSize = label_font_size;

% Loop through each of the experimental sessions
for xx = 1:length(avg_depth)

    % Scatter of depth of modulation
    scatter(target_centers{xx}, avg_depth{xx}, marker_size, marker_shape, 'LineWidth', 1.5);

    % Line of best fit
    if isequal(error_choice, 1)
        best_fit = fitlm(target_centers{xx}, avg_depth{xx});
    
        % Slopes and y_intercepts
        units_slope_and_y_intercept(1) = table2array(best_fit.Coefficients(2,1));
        units_slope_and_y_intercept(2) = table2array(best_fit.Coefficients(1,1));
    
        best_fit = refline(units_slope_and_y_intercept(1), units_slope_and_y_intercept(2));
        best_fit.LineWidth = 1;
    end

end % End the xds loop

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = numel(findobj('type','figure')):-1:1
        title_string(ii) = strrep(title_string(ii), ':', '');
        title_string(ii) = strrep(title_string(ii), '.0', '');
        title_string(ii) = strrep(title_string(ii), 'vs.', 'vs');
        title_string(ii) = strrep(title_string(ii), 'mg.', 'mg');
        title_string(ii) = strrep(title_string(ii), 'kg.', 'kg');
        title_string(ii) = strrep(title_string(ii), '.', '_');
        title_string(ii) = strrep(title_string(ii), '/', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(title_string(ii))), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(title_string(ii))), 'png')
            saveas(gcf, fullfile(save_dir, char(title_string(ii))), 'pdf')
            saveas(gcf, fullfile(save_dir, char(title_string(ii))), 'fig')
        end
        close gcf
    end
end







