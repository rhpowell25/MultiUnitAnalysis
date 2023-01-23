
%% Load the output structures
clc

% What Drug Do You Want To Load? ('Caff', 'Cyp', 'Lex', 'Con') or ('YYYYMMDD')
Drug_Choice = 'Caff';
% What Monkey Do You Want To Load? ('All', 'Mihili', 'Jango', 'Jaco', 'Pop', 'Pancake')
Monkey_Choice = strings;
Monkey_Choice{1} = 'Pop';
%Monkey_Choice{1} = 'Pancake';
% What targets do you want to use ('Max', 'Min')
Tgt_Choice = 'Max';

event = 'window_trial_gocue_2';

Sampling_Params = struct( ...
    'unit_quality', 'Stable', ... % Unit quality ('All' vs. 'Stable')
    'ISI_quality', 'All', ... % ISI quality ('All' vs. 'Single')
    'depth_change', NaN, ... % What change in depth of modulation do you want to observe (# vs NaN)
    'pref_dir', 'All', ... % What preferred direction do you want to plot(-90, 0, 90, 180, 'All')
    'depth_min', NaN, ... % What minimum depth of modulation do you want to observe (# vs NaN)
    'depth_sig', 'All'); % Unit modulation significance ('All' vs. 'Sig')

[bsfr_morn, bsfr_noon, depth_morn, depth_noon, ~, ~, ~, ~, ~, ~, wave_p_value, nonlin_p_value, fract_contam, ...
    pref_dir, num_targets, drug_dose, file_name, all_unit_names] = Fast_Load_Depth(Drug_Choice, Monkey_Choice, event, Tgt_Choice, Sampling_Params);

%% Some of the plotting specifications
% Individual Sessions or All Sessions? ('per_trial' vs 'all_trials')
trial_sessions = 'all_trials';

% Plot the depth or the baseline ('Depth' vs. 'Baseline' vs 'Both'?
depth_vs_bsfr = 'Baseline';

% Which task do you want to plot ('PG' vs. 'WS', vs ''All')?
trial_task = 'Both';

% Do you want to plot the morning / afternoon legend? (1 = Yes, 0 = No)
plot_legend = 1;

% Remove sessions with only a single target? (1 = Yes, 0 = No)
max_vs_min = 0;

% What statistical measure do you want to use ('T-Test', 'Wilcox')
stat_test = 'T-Test';

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
Save_Figs = 0;
if ~isequal(Save_Figs, 0)
    % Do you want a save title or blank title (1 = save_title, 0 = blank)
    Fig_Save_Title = 0;
    close all
end

%% Some variable extraction & definitions

% Font specifications
label_font_size = 30;
legend_size = 30;
mean_line_width = 5;
if isequal(plot_legend, 1)
    p_value_dims = [0.51 0.3 0.44 0.44];
else
    p_value_dims = [0.51 0.45 0.44 0.44];
end
if isequal(plot_legend, 1)
    n_value_dims = [0.49 0.225 0.44 0.44];
else
    n_value_dims = [0.49 0.375 0.44 0.44];
end
title_font_size = 14;
save_title_font_size = 45;
font_name = 'Arial';

% Save Counter
ss = 1;

if ~isequal(Save_Figs, 0)
    save_title = strings;
end

%% Remove experiments with only a single target per direction
if isequal(max_vs_min, 1)
    disp('Removing single target sessions')
    for ii = 1:length(file_name)
        single_tgt_idx = find(ismember(num_targets{ii}, 1));
        all_unit_names{ii}(single_tgt_idx) = [];
        bsfr_morn{ii}(single_tgt_idx) = [];
        depth_morn{ii}(single_tgt_idx) = [];
        bsfr_noon{ii}(single_tgt_idx) = [];
        depth_noon{ii}(single_tgt_idx) = [];
        fract_contam{ii}(single_tgt_idx) = [];
        pref_dir{ii}(single_tgt_idx) = [];
    end
end

%% Only use trial tasks selected
if strcmp(trial_task, 'PG')
    disp('Powergrasp Only')
    powergrasp_idx = find(~contains(file_name, 'PG'));
    all_unit_names(powergrasp_idx) = [];
    bsfr_morn(powergrasp_idx) = [];
    depth_morn(powergrasp_idx) = [];
    bsfr_noon(powergrasp_idx) = [];
    depth_noon(powergrasp_idx) = [];
    fract_contam(powergrasp_idx) = [];
    pref_dir(powergrasp_idx) = [];
    drug_dose(powergrasp_idx) = [];
    file_name(powergrasp_idx) = [];
end

if strcmp(trial_task, 'WS')
    disp('Wrist w/ Spring Only')
    wrist_idx = find(~contains(file_name, 'WS'));
    all_unit_names(wrist_idx) = [];
    bsfr_morn(wrist_idx) = [];
    depth_morn(wrist_idx) = [];
    bsfr_noon(wrist_idx) = [];
    depth_noon(wrist_idx) = [];
    fract_contam(wrist_idx) = [];
    pref_dir(wrist_idx) = [];
    drug_dose(wrist_idx) = [];
    file_name(wrist_idx) = [];
end

%% Merge all the information (Only if you selected 'all_trial')
if strcmp(trial_sessions, 'all_trials')
    % Find the amount of units each experiment uses
    units_per_experiment = zeros(length(all_unit_names),1);
    for ii = 1:length(units_per_experiment)
        units_per_experiment(ii) = length(all_unit_names{ii});
    end
    total_units = sum(units_per_experiment);

    % Concatenate all the information
    all_trials_unit_names = struct([]);
    all_trials_bsfr_morn = zeros(length(total_units),1);
    all_trials_depth_morn = zeros(length(total_units),1);
    all_trials_bsfr_noon = zeros(length(total_units),1);
    all_trials_depth_noon = zeros(length(total_units),1);
    all_trials_merged_fract_contam = zeros(length(total_units),1);
    cc = 1;
    for xx = 1:length(all_unit_names)
        for jj = 1:length(all_unit_names{xx})
            all_trials_unit_names{cc,1} = all_unit_names{xx}(jj);
            all_trials_bsfr_morn(cc,1) = bsfr_morn{xx}(jj);
            all_trials_depth_morn(cc,1) = depth_morn{xx}(jj);
            all_trials_bsfr_noon(cc,1) = bsfr_noon{xx}(jj);
            all_trials_depth_noon(cc,1) = depth_noon{xx}(jj);
            all_trials_merged_fract_contam(cc,1) = fract_contam{xx}(jj);
            cc = cc + 1;
        end
    end
    % Rename the now merged variables
    all_unit_names = struct([]);
    all_unit_names{1,1} = all_trials_unit_names;
    bsfr_morn = struct([]);
    bsfr_morn{1,1} = all_trials_bsfr_morn;
    depth_morn = struct([]);
    depth_morn{1,1} = all_trials_depth_morn;
    bsfr_noon = struct([]);
    bsfr_noon{1,1} = all_trials_bsfr_noon;
    depth_noon = struct([]);
    depth_noon{1,1} = all_trials_depth_noon;
    fract_contam = struct([]);
    fract_contam{1,1} = all_trials_merged_fract_contam;
end

%% Merge all the information (Only if you selected 'all_trial')
if strcmp(trial_sessions, 'all_trials')
    % Find the amount of units each experiment uses
    units_per_experiment = zeros(length(all_unit_names),1);
    for ii = 1:length(units_per_experiment)
        units_per_experiment(ii) = length(all_unit_names{ii});
    end
    total_units = sum(units_per_experiment);

    % Concatenate all the information
    all_trials_unit_names = struct([]);
    all_trials_bsfr_morn = zeros(length(total_units),1);
    all_trials_depth_morn = zeros(length(total_units),1);
    all_trials_bsfr_noon = zeros(length(total_units),1);
    all_trials_depth_noon = zeros(length(total_units),1);
    all_trials_merged_fract_contam = zeros(length(total_units),1);
    cc = 1;
    for xx = 1:length(all_unit_names)
        for jj = 1:length(all_unit_names{xx})
            all_trials_unit_names{cc,1} = all_unit_names{xx}(jj);
            all_trials_bsfr_morn(cc,1) = bsfr_morn{xx}(jj);
            all_trials_depth_morn(cc,1) = depth_morn{xx}(jj);
            all_trials_bsfr_noon(cc,1) = bsfr_noon{xx}(jj);
            all_trials_depth_noon(cc,1) = depth_noon{xx}(jj);
            all_trials_merged_fract_contam(cc,1) = fract_contam{xx}(jj);
            cc = cc + 1;
        end
    end
    % Rename the now merged variables
    all_unit_names = struct([]);
    all_unit_names{1,1} = all_trials_unit_names;
    bsfr_morn = struct([]);
    bsfr_morn{1,1} = all_trials_bsfr_morn;
    depth_morn = struct([]);
    depth_morn{1,1} = all_trials_depth_morn;
    bsfr_noon = struct([]);
    bsfr_noon{1,1} = all_trials_bsfr_noon;
    depth_noon = struct([]);
    depth_noon{1,1} = all_trials_depth_noon;
    fract_contam = struct([]);
    fract_contam{1,1} = all_trials_merged_fract_contam;
end

%% Run through each experiment seperately

for xx = 1:length(all_unit_names)

    % Skip the function if no units match
    if isempty(all_unit_names{xx})
        disp('No units match criteria')
        continue
    end

    %% Add the monkey name to the title

    if strcmp(Monkey_Choice, 'All')
        fig_title = strcat('All Monkeys,', {' '});
    else
        fig_title = strcat(Monkey_Choice, {' '});
    end

    %% Depth of modulation vs. Baseline firing rate
    if strcmp(depth_vs_bsfr, 'Both')
        depth_base = 2;
    else
        depth_base = 1;
    end

    %% Reassign variables according to what you're plotting

    for db = 1:depth_base

        if strcmp(depth_vs_bsfr, 'Baseline')
            disp('Baseline Firing Rate')
            fire_rate_morn = bsfr_morn;
            fire_rate_noon = bsfr_noon;
        end
        if strcmp(depth_vs_bsfr, 'Depth')
            disp('Depth of Modulation')
            fire_rate_morn = depth_morn;
            fire_rate_noon = depth_noon;
        end

        if strcmp(depth_vs_bsfr, 'Both')
            if isequal(db, 1)
                disp('Baseline Firing Rate')
                fire_rate_morn = bsfr_morn;
                fire_rate_noon = bsfr_noon;
            end
            if isequal(db, 2)
                disp('Depth of Modulation')
                fire_rate_morn = depth_morn;
                fire_rate_noon = depth_noon;
            end
        end

        %% Do the statistics

        if strcmp(stat_test, 'T-Test')
            [~, fire_rate_p_val] = ttest(fire_rate_morn{xx,1}, fire_rate_noon{xx,1});
        elseif strcmp(stat_test, 'Wilcox')
            [fire_rate_p_val, ~] = ranksum(fire_rate_morn{xx,1}, fire_rate_noon{xx,1});
        end
        fire_rate_n_val = length(fire_rate_morn{xx,1});

        %% Plot the depth of modulation histograms

        % Find the mean baseline firing rates
        fire_rate_mean_morn = mean(fire_rate_morn{xx,1});
        fire_rate_mean_noon = mean(fire_rate_noon{xx,1});
        fire_rate_perc_change = ((fire_rate_mean_noon - fire_rate_mean_morn) / fire_rate_mean_morn) * 100;

        % Merged
        fire_rate_fig = figure;
        fire_rate_fig.Position = [200 50 800 600];
        hold on

        histogram(fire_rate_morn{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
        histogram(fire_rate_noon{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])

        fprintf('The effect size is %0.2f percent \n', fire_rate_perc_change)

        % Set the axis
        x_limits = xlim;
        y_limits = ylim;
        xlim([x_limits(1), x_limits(2) + 5])
        ylim([y_limits(1), y_limits(2) + 0.25])

        % Plot dummy points for the legend
        if isequal(plot_legend, 1)
            dummy_morn = plot(-1,-1, 's', 'MarkerSize',20, 'LineWidth', 1.5, ...
                    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], 'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
            dummy_noon = plot(-2,-1, 's', 'MarkerSize',20, 'LineWidth', 1.5, ...
                    'MarkerEdgeColor',[.5 0 .5], 'MarkerFaceColor',[.5 0 .5]);
        end

        % Plot the means
        line([fire_rate_mean_morn fire_rate_mean_morn], [y_limits(1) y_limits(2) + 0.25], ... 
            'LineStyle','--', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', mean_line_width)
        line([fire_rate_mean_noon fire_rate_mean_noon], [y_limits(1) y_limits(2) + 0.25], ... 
            'LineStyle','--', 'Color', [.5 0 .5], 'LineWidth', mean_line_width)
            
        % Set the title
        if strcmp(trial_sessions, 'per_trial')
            if ~strcmp(Drug_Choice, 'Con')
                title(strcat(fig_title, file_name{xx}, {' '}, drug_dose(xx)), 'FontSize', title_font_size)
            else
                title(strcat(fig_title, file_name{xx}), 'FontSize', title_font_size)
            end
        elseif strcmp(trial_sessions, 'all_trials')
            scatter_title = strcat('All Trials,', {' '}, Drug_Choice);
            if strcmp(trial_task, 'PG')
                scatter_title = strcat(scatter_title, {' '}, 'PG');
            end
            if strcmp(trial_task, 'WS')
                scatter_title = strcat(scatter_title, {' '}, 'WS');
            end
            title(strcat(fig_title, scatter_title, {' '}, 'Depth'), 'FontSize', title_font_size)
        end

        % Axis Editing
        figure_axes = gca;
        % Set ticks to outside
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set the tick label font size
        figure_axes.FontSize = label_font_size - 15;

        % Label the axis
        ylabel('Units', 'FontSize', label_font_size)
        if strcmp(depth_vs_bsfr, 'Baseline')
            xlabel('Baseline Firing Rate (Hz)', 'FontSize', label_font_size)
        end
        if strcmp(depth_vs_bsfr, 'Depth')
            xlabel('Depth of Modulation (Hz)', 'FontSize', label_font_size)
        end
        if strcmp(depth_vs_bsfr, 'Both')
            if isequal(db, 1)
                xlabel('Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
            end
            if isequal(db, 2)
                xlabel('Depth of Modulation (Hz)', 'FontSize', label_font_size);
            end
        end

        if ~isequal(Save_Figs, 0)
            fig_info = get(gca,'title');
            save_title(ss) = get(fig_info, 'string');
            % Make the title the drug
            if strcmp(Drug_Choice, 'Lex')
                title_color = [0 0.5 0];
                fig_title = 'Escitalopram';
            end
            if strcmp(Drug_Choice, 'Caff')
                title_color = [0 0.5 0];
                fig_title = 'Caffeine';
            end
            if strcmp(Drug_Choice, 'Cyp')
                title_color =  'r';
                fig_title = 'Cyproheptadine';
            end
            if strcmp(Drug_Choice, 'Con')
                title_color =  'k';
                fig_title = 'Control';
            end
            if contains(Drug_Choice, '202')
                title_color =  'k';
                fig_title = strcat(fig_title, file_name{xx});
            end
            if ~isequal(Fig_Save_Title, 0)
                if strcmp(trial_task, 'All')
                    title(strcat(fig_title, fig_title), 'FontSize', save_title_font_size - 15)
                else
                    save_title(ss) = get(fig_info, 'string');
                end
                %title(fig_title, 'Fontsize', save_title_font_size, 'Color', title_color);
            else
                title('')
            end
        end

        % Annotation of the p_value
        if round(fire_rate_p_val, 3) > 0
            p_value_string = strcat('p =', {' '}, mat2str(round(fire_rate_p_val, 3)));
            p_value_string = {char(p_value_string)};
            ann_p_value = annotation('textbox', p_value_dims, 'String', p_value_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_p_value.FontSize = legend_size;
            ann_p_value.FontName = font_name;
        end

        if isequal(round(fire_rate_p_val, 3), 0)
            p_value_string = strcat('p <', {' '}, '0.001');
            p_value_string = {char(p_value_string)};
            ann_p_value = annotation('textbox', p_value_dims, 'String', p_value_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_p_value.FontSize = legend_size;
            ann_p_value.FontName = font_name;
        end

        % Annotation of the n_value
        n_value_string = strcat('n =', {' '}, mat2str(round(fire_rate_n_val, 3)));
        n_value_string = {char(n_value_string)};
        ann_n_value = annotation('textbox', n_value_dims, 'String', n_value_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_n_value.FontSize = legend_size;
        ann_n_value.FontName = font_name;

        if isequal(plot_legend, 1)
            legend([dummy_morn, dummy_noon], ... 
                {'Morning', 'Afternoon'}, ... 
                'FontSize', legend_size, 'Location', 'NorthEast')
            legend boxoff
        end

        % Only label every other tick
        %x_labels = string(figure_axes.XAxis.TickLabels);
        %y_labels = string(figure_axes.YAxis.TickLabels);
        %x_labels(2:2:end) = NaN;
        %y_labels(2:2:end) = NaN;
        %figure_axes.XAxis.TickLabels = x_labels;
        %figure_axes.YAxis.TickLabels = y_labels;

        % Add to the counter
        ss = ss + 1;

    end
end

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = numel(findobj('type','figure')):-1:1
        save_title(ii) = strrep(save_title(ii), ':', '');
        save_title(ii) = strrep(save_title(ii), '.0', '');
        save_title(ii) = strrep(save_title(ii), 'vs.', 'vs');
        save_title(ii) = strrep(save_title(ii), 'mg.', 'mg');
        save_title(ii) = strrep(save_title(ii), 'kg.', 'kg');
        save_title(ii) = strrep(save_title(ii), '.', '_');
        save_title(ii) = strrep(save_title(ii), '/', '_');
        save_title(ii) = strcat(save_title(ii), {' '}, '(Hist)');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), 'png')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), 'pdf')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), 'fig')
        end
        close gcf
    end
end


