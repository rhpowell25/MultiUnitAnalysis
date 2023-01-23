
%% Load the output structures
% What Drug Do You Want To Load? ('Caff', 'Cyp', 'Lex', 'Con') or ('YYYYMMDD')
Drug_Choice = 'Cyp';
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
% Individual Sessions or All Sessions? ('per_trial' vs 'all_trials'
trial_sessions = 'per_trial';

% Plot the depth or the baseline ('Depth' vs. 'Baseline' vs 'Both'?
depth_vs_bsfr = 'Depth';

% Which task do you want to plot ('PG' vs. 'WS', vs 'All')?
trial_task = 'All';

% What range of depth of modulation do you want to exclude
depth_exclusion = [-5,5];

% Do you want to remove outliers? (1 = Yes, 0 = No)
remove_outliers = 1;

% Do you want to manually set the y-axis?
man_y_axis = 'No';
%man_y_axis = [-300, 300];

% Remove sessions with only a single target?(1 = Yes, 0 = No)
max_vs_min = 0;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
Save_Figs = 0;
if ~isequal(Save_Figs, 0)
    close all
end

%% Some variable extraction & definitions

% Font specifications
label_font_size = 30;
zero_line_width = 2;
legend_size = 20;
effect_size_dims = [0.7 0.375 0.44 0.44];
if isequal(remove_outliers, 1)
    outlier_ann_dims = [0.7 0.325 0.44 0.44];
end
title_font_size = 14;
save_title_font_size = 45;
font_name = 'Arial';
fig_size = 700;

% Save Counter
ss = 1;

if ~isequal(Save_Figs, 0)
    save_title = strings;
    % Do you want a save title or blank title (1 = save_title, 0 = blank)
    Fig_Save_Title = 0;
    close all
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
    all_trials_unit_names = strings;
    all_trials_bsfr_morn = zeros(length(total_units),1);
    all_trials_depth_morn = zeros(length(total_units),1);
    all_trials_bsfr_noon = zeros(length(total_units),1);
    all_trials_depth_noon = zeros(length(total_units),1);
    all_trials_merged_fract_contam = zeros(length(total_units),1);
    cc = 1;
    for xx = 1:length(all_unit_names)
        for jj = 1:length(all_unit_names{xx})
            all_trials_unit_names(cc,1) = all_unit_names{xx}(jj);
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

%% Remove any units that fall below the minimum selected depth of modulation

if ~isnan(depth_exclusion)
    depth_violations_morn = struct([]);
    depth_violations_noon = struct([]);
    depth_violations = struct([]);
    for ii = 1:length(all_unit_names)
        depth_violations_morn{ii,1} = find(depth_morn{ii,1} > depth_exclusion(1) & depth_morn{ii,1} < depth_exclusion(2));
        depth_violations_noon{ii,1} = find(depth_noon{ii,1} > depth_exclusion(1) & depth_noon{ii,1} < depth_exclusion(2));
        depth_violations{ii,1} = intersect(depth_violations_morn{ii,1}, depth_violations_noon{ii,1});
        all_unit_names{ii,1}(depth_violations{ii,1}) = [];
        bsfr_morn{ii,1}(depth_violations{ii,1}) = [];
        depth_morn{ii,1}(depth_violations{ii,1}) = [];
        bsfr_noon{ii,1}(depth_violations{ii,1}) = [];
        depth_noon{ii,1}(depth_violations{ii,1}) = [];
        fract_contam{ii,1}(depth_violations{ii,1}) = [];
    end
end

%% Loop through each of the experimental sessions
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

        %% Plot the Violin Plot

        effect_sizes = (fire_rate_noon{xx,1} - fire_rate_morn{xx,1}) ./ abs(fire_rate_morn{xx,1}) * 100;

        % Remove outliers
        if isequal(remove_outliers, 1)
            effect_outliers = isoutlier(effect_sizes, 'mean');
            %if ~isnan(depth_exclusion)
            %    num_outliers = length(find(effect_outliers == 1)) + length(depth_violations{xx,1});
            %elseif isnan(depth_exclusion)
                num_outliers = length(find(effect_outliers == 1));
            %end
            effect_sizes(effect_outliers) = [];
        end

        if strcmp(Drug_Choice, 'Caff') || strcmp(Drug_Choice, 'Lex')
            violin_color = [0 0.5 0];
        elseif strcmp(Drug_Choice, 'Cyp')
            violin_color = [1 0 0];
        else
            violin_color = [0 0 0];
        end

        violin_fig = figure;
        violin_fig.Position = [200 50 fig_size fig_size];
        hold on
        effect_positions = (1:length(effect_sizes));
        Violin_Plot({effect_sizes}, effect_positions, 'ViolinColor', violin_color);

        % Set the axis
        x_limits = xlim;
        y_limits = ylim;
        if ~ischar(man_y_axis)
            ylim(man_y_axis)
        end
        xlim([0.5, 1.5])

        ylabel('Effect Size (%)', 'FontSize', label_font_size)

        % Draw the unity line 
        line([0.5, 1.5], [0, 0], 'Color', 'k', 'Linewidth', zero_line_width, 'Linestyle','--')

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
        figure_axes.FontSize = label_font_size - 10;

        if ~isequal(Save_Figs, 0)
            fig_info = get(gca,'title');
            save_title(ss) = get(fig_info, 'string');
            % Make the title the drug
            if strcmp(Drug_Choice, 'Lex')
                fig_title = 'Escitalopram';
            end
            if strcmp(Drug_Choice, 'Caff')
                fig_title = 'Caffeine';
            end
            if strcmp(Drug_Choice, 'Cyp')
                fig_title = 'Cyproheptadine';
            end
            if strcmp(Drug_Choice, 'Con')
                fig_title = 'Control';
            end
            if contains(Drug_Choice, '202')
                fig_title = strcat(units_title, file_name{xx});
            end
            if ~isequal(Fig_Save_Title, 0)
                if strcmp(trial_task, 'All')
                    title(strcat(units_title, fig_title), 'FontSize', save_title_font_size - 15)
                else
                    save_title(ss) = get(fig_info, 'string');
                end
                title(fig_title, 'Fontsize', save_title_font_size, 'Color', violin_color);
            else
                title('')
            end
        end

        % Only label every other tick
        x_labels = string(figure_axes.XAxis.TickLabels);
        y_labels = string(figure_axes.YAxis.TickLabels);
        x_labels(1:end) = NaN;
        y_labels(2:2:end) = NaN;
        figure_axes.XAxis.TickLabels = x_labels;
        figure_axes.YAxis.TickLabels = y_labels;

        % Display the percent change & statistics
        fire_rate_mean_morn = mean(fire_rate_morn{xx,1});
        fire_rate_mean_noon = mean(fire_rate_noon{xx,1});
        [~, fire_rate_p_val] = ttest(fire_rate_morn{xx,1}, fire_rate_noon{xx,1});
        fire_rate_perc_change = ((fire_rate_mean_noon - fire_rate_mean_morn) / abs(fire_rate_mean_morn)) * 100;
        fprintf('Percent change is %0.2f%%, p = %0.3f \n', fire_rate_perc_change, fire_rate_p_val)

        % Annotation of the effect size
        effect_size_string = strcat('Δ% =', {' '}, mat2str(round(fire_rate_perc_change, 1)), '%');
        effect_ann_string = {char(effect_size_string)};
        effect_ann = annotation('textbox', effect_size_dims, 'String', effect_ann_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'Left');
        effect_ann.FontSize = legend_size;
        effect_ann.FontName = font_name;

        % Annotation of the outliers removed
        if isequal(remove_outliers, 1)
            if ~ischar(man_y_axis)
                num_outliers = num_outliers + length(find(effect_sizes > man_y_axis(2))) + ...
                    length(find(effect_sizes < man_y_axis(1)));
            end
            outlier_string = strcat(mat2str(num_outliers), {' '}, '> 3σ');
            outlier_ann_string = {char(outlier_string)};
            outlier_ann = annotation('textbox', outlier_ann_dims, 'String', outlier_ann_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'Left');
            outlier_ann.FontSize = legend_size;
            outlier_ann.FontName = font_name;
        end


        % Add to the counter
        ss = ss + 1;

    end

end % End the xds loop

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
        save_title(ii) = strcat(save_title(ii), {' '}, '(Violin)');
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







