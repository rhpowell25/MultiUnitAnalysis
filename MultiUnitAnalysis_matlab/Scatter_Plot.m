
%% Load the output structures
% What Drug Do You Want To Load? ('Caff', 'Cyp', 'Lex', 'Con') or ('YYYYMMDD')
Drug_Choice = 'Con';
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

[bsfr_morn, bsfr_noon, depth_morn, depth_noon, bsfr_err_morn, mpfr_err_morn, bsfr_err_noon, mpfr_err_noon, ...
    depth_t_test, depth_wilcoxon, wave_p_value, nonlin_p_value, fract_contam, pref_dir, num_targets, ...
    drug_dose, file_name, all_unit_names] = Fast_Load_Depth(Drug_Choice, Monkey_Choice, event, Tgt_Choice, Sampling_Params);

%% Some of the plotting specifications
% Individual Sessions or All Sessions? ('per_trial' vs 'all_trials'
trial_sessions = 'per_trial';

% Plot the depth or the baseline ('Depth' vs. 'Baseline' vs 'Both'?
depth_vs_bsfr = 'Depth';

% Which task do you want to plot ('PG' vs. 'WS', vs 'All')?
trial_task = 'PG';

% Best fit or elliptial error probable? ('none', 'best_fit', or 'ellip_err_prob')
error_choice = 'none';

% Which stability metric do you want to use ('Wave', "Nonlin')
stability_choice = 'Nonlin';

% Which significance metric do you want to use ('T-Test', 'Wilcox')
sig_choice = 'T-Test';

% Do you want the name of each unit labeled? (1 = Yes, 0 = No)
unit_label = 0;

% Do you want to plot the hypothesis arrows (1 = Yes, 0 = No)
hypo_arrows = 0;

% Remove sessions with only a single target?(1 = Yes, 0 = No)
max_vs_min = 0;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
Save_Figs = 'png';
if ~isequal(Save_Figs, 0)
    close all
end

%% Some variable extraction & definitions

% Font specifications
label_font_size = 30;
legend_size = 30;
title_font_size = 14;
save_title_font_size = 45;
font_name = 'Arial';
fig_size = 700;

% Save Counter
ss = 1;

if ~isequal(Save_Figs, 0)
    save_title = strings;
end

% Scatter Marker Shapes
single_marker = '.';
multi_marker = '*';

% Scatter Marker Colors
insig_color = 0;
sig_color = 1;

% Scatter Marker sizes
single_marker_size = 1000;
multi_marker_size = 250;

%% Use chosen stability metric

if strcmp(stability_choice, 'Nonlin')
    stable_metric = nonlin_p_value;
elseif strcmp(stability_choice, 'Wave')
    stable_metric = wave_p_value;
end

%% Use chosen significance metric

if strcmp(sig_choice, 'T-Test')
    sig_metric = depth_t_test;
elseif strcmp(sig_choice, 'Wilcox')
    sig_metric = depth_wilcoxon;
end

%% Remove experiments with only a single target per direction
if isequal(max_vs_min, 1)
    disp('Removing single target sessions')
    for ii = 1:length(file_name)
        single_tgt_idx = find(ismember(num_targets{ii}, 1));
        all_unit_names{ii}(single_tgt_idx) = [];
        bsfr_morn{ii}(single_tgt_idx) = [];
        bsfr_err_morn{ii}(single_tgt_idx) = [];
        depth_morn{ii}(single_tgt_idx) = [];
        mpfr_err_morn{ii}(single_tgt_idx) = [];
        bsfr_noon{ii}(single_tgt_idx) = [];
        bsfr_err_noon{ii}(single_tgt_idx) = [];
        depth_noon{ii}(single_tgt_idx) = [];
        mpfr_err_noon{ii}(single_tgt_idx) = [];
        stable_metric{ii}(single_tgt_idx) = [];
        sig_metric{ii}(single_tgt_idx) = [];
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
    bsfr_err_morn(powergrasp_idx) = [];
    depth_morn(powergrasp_idx) = [];
    mpfr_err_morn(powergrasp_idx) = [];
    bsfr_noon(powergrasp_idx) = [];
    bsfr_err_noon(powergrasp_idx) = [];
    depth_noon(powergrasp_idx) = [];
    mpfr_err_noon(powergrasp_idx) = [];
    stable_metric(powergrasp_idx) = [];
    sig_metric(powergrasp_idx) = [];
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
    bsfr_err_morn(wrist_idx) = [];
    depth_morn(wrist_idx) = [];
    mpfr_err_morn(wrist_idx) = [];
    bsfr_noon(wrist_idx) = [];
    bsfr_err_noon(wrist_idx) = [];
    depth_noon(wrist_idx) = [];
    mpfr_err_noon(wrist_idx) = [];
    stable_metric(wrist_idx) = [];
    sig_metric(wrist_idx) = [];
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
    all_trials_bsfr_err_morn = zeros(length(total_units),1);
    all_trials_depth_morn = zeros(length(total_units),1);
    all_trials_mpfr_err_morn = zeros(length(total_units),1);
    all_trials_bsfr_noon = zeros(length(total_units),1);
    all_trials_bsfr_err_noon = zeros(length(total_units),1);
    all_trials_depth_noon = zeros(length(total_units),1);
    all_trials_mpfr_err_noon = zeros(length(total_units),1);
    all_trials_sort_p_value = zeros(length(total_units),1);
    all_trials_sig_p_value = zeros(length(total_units),1);
    all_trials_merged_fract_contam = zeros(length(total_units),1);
    cc = 1;
    for xx = 1:length(all_unit_names)
        for jj = 1:length(all_unit_names{xx})
            all_trials_unit_names(cc,1) = all_unit_names{xx}(jj);
            all_trials_bsfr_morn(cc,1) = bsfr_morn{xx}(jj);
            all_trials_bsfr_err_morn(cc,1) = bsfr_err_morn{xx}(jj);
            all_trials_depth_morn(cc,1) = depth_morn{xx}(jj);
            all_trials_mpfr_err_morn(cc,1) = mpfr_err_morn{xx}(jj);
            all_trials_bsfr_noon(cc,1) = bsfr_noon{xx}(jj);
            all_trials_bsfr_err_noon(cc,1) = bsfr_err_noon{xx}(jj);
            all_trials_depth_noon(cc,1) = depth_noon{xx}(jj);
            all_trials_mpfr_err_noon(cc,1) = mpfr_err_noon{xx}(jj);
            all_trials_sort_p_value(cc,1) = stable_metric{xx}(jj);
            all_trials_sig_p_value(cc,1) = sig_metric{xx}(jj);
            all_trials_merged_fract_contam(cc,1) = fract_contam{xx}(jj);
            cc = cc + 1;
        end
    end
    % Rename the now merged variables
    all_unit_names = struct([]);
    all_unit_names{1,1} = all_trials_unit_names;
    bsfr_morn = struct([]);
    bsfr_morn{1,1} = all_trials_bsfr_morn;
    bsfr_err_morn = struct([]);
    bsfr_err_morn{1,1} = all_trials_bsfr_err_morn;
    depth_morn = struct([]);
    depth_morn{1,1} = all_trials_depth_morn;
    mpfr_err_morn = struct([]);
    mpfr_err_morn{1,1} = all_trials_mpfr_err_morn;
    bsfr_noon = struct([]);
    bsfr_noon{1,1} = all_trials_bsfr_noon;
    bsfr_err_noon = struct([]);
    bsfr_err_noon{1,1} = all_trials_bsfr_err_noon;
    depth_noon = struct([]);
    depth_noon{1,1} = all_trials_depth_noon;
    mpfr_err_noon = struct([]);
    mpfr_err_noon{1,1} = all_trials_mpfr_err_noon;
    stable_metric = struct([]);
    stable_metric{1,1} = all_trials_sort_p_value;
    sig_metric = struct([]);
    sig_metric{1,1} = all_trials_sig_p_value;
    fract_contam = struct([]);
    fract_contam{1,1} = all_trials_merged_fract_contam;
end

%% Loop through each of the experimental sessions
for xx = 1:length(all_unit_names)

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
            fire_rate_err_morn = bsfr_err_morn;
            fire_rate_noon = bsfr_noon;
            fire_rate_err_noon = bsfr_err_noon;
        end
        if strcmp(depth_vs_bsfr, 'Depth')
            disp('Depth of Modulation')
            fire_rate_morn = depth_morn;
            fire_rate_err_morn = mpfr_err_morn;
            fire_rate_noon = depth_noon;
            fire_rate_err_noon = mpfr_err_noon;
        end

        if strcmp(depth_vs_bsfr, 'Both')
            if isequal(db, 1)
                disp('Baseline Firing Rate')
                fire_rate_morn = bsfr_morn;
                fire_rate_err_morn = bsfr_err_morn;
                fire_rate_noon = bsfr_noon;
                fire_rate_err_noon = bsfr_err_noon;
            end
            if isequal(db, 2)
                disp('Depth of Modulation')
                fire_rate_morn = depth_morn;
                fire_rate_err_morn = mpfr_err_morn;
                fire_rate_noon = depth_noon;
                fire_rate_err_noon = mpfr_err_noon;
            end
        end

        %% Extract the lines of best fit

        if strcmp(error_choice, 'best_fit')
            % Line of best fit
            units_best_fit = fitlm(fire_rate_morn{xx}, fire_rate_noon{xx});

            % Slopes and y_intercepts
            units_slope_and_y_intercept(1) = table2array(units_best_fit.Coefficients(2,1));
            units_slope_and_y_intercept(2) = table2array(units_best_fit.Coefficients(1,1));

            % Find the r^2 & square standard dev
            units_r_squared = units_best_fit.Rsquared.Ordinary;
            units_least_square_std = sqrt(units_best_fit.SSE / (length(fire_rate_morn{xx,1}) - 2));
        elseif strcmp(error_choice, 'ellip_err_prob') && length(fire_rate_morn{xx}) > 2
            % Plot the elliptical error probable
            scatter_x_axis = fire_rate_morn{xx}';
            scatter_y_axis = fire_rate_noon{xx}';
            % The percent of points you want
            err_percent = .5;
            % Run the elliptical error probability function
            [~, ~, ~, ~, ~, X_ellipse, Y_ellipse] = ... 
                Ellip_Err_Prob(scatter_x_axis, scatter_y_axis, err_percent);
        end

        %% Plot the Depth of Modulation Scatter

        scatter_fig = figure;
        scatter_fig.Position = [200 50 fig_size fig_size];
        hold on

        % Plot a circle (x-a) = rcos(t) / (y-b) = rsin(t)
        %t = linspace(0,2*pi);
        %rad = 7;
        %x = rad*cos(t) + fire_rate_morn{1,1}(12);
        %y = rad*sin(t) + fire_rate_noon{1,1}(12);
        %plot(x,y, 'r', 'LineWidth', 2)

        % Set the title
        if strcmp(trial_sessions, 'per_trial')
            title(strcat(fig_title, file_name{xx}, {' '}, drug_dose(xx)), 'FontSize', title_font_size)
        elseif strcmp(trial_sessions, 'all_trials')
            scatter_title = 'All Trials';
            title(strcat(fig_title, scatter_title, {' '}, Drug_Choice), 'FontSize', title_font_size)
        end

        % Axis Editing
        figure_axes = gca;
        % Set ticks to outside
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set the tick label font size
        figure_axes.FontSize = label_font_size - 5;

        % Label the axis
        if strcmp(depth_vs_bsfr, 'Baseline')
            xlabel('Morning Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
            ylabel('Afternoon Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
        end
        if strcmp(depth_vs_bsfr, 'Depth')
            xlabel('Morning Depth of Modulation (Hz)', 'FontSize', label_font_size);
            ylabel('Afternoon Depth of Modulation (Hz)', 'FontSize', label_font_size);
        end

        if strcmp(depth_vs_bsfr, 'Both')
            if isequal(db, 1)
                xlabel('Morning Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
                ylabel('Afternoon Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
            end
            if isequal(db, 2)
                xlabel('Morning Depth of Modulation (Hz)', 'FontSize', label_font_size);
                ylabel('Afternoon Depth of Modulation (Hz)', 'FontSize', label_font_size);
            end
        end

        % Calculate the axis limits
        min_err_morn = min(fire_rate_morn{xx} - fire_rate_err_morn{xx});
        min_err_noon = min(fire_rate_noon{xx} - fire_rate_err_noon{xx});
        axis_min = round(min(min_err_morn, min_err_noon)/5)*5;
        max_err_morn = max(fire_rate_morn{xx} + fire_rate_err_morn{xx});
        max_err_noon = max(fire_rate_noon{xx} + fire_rate_err_noon{xx});
        axis_max = round(max(max_err_morn, max_err_noon)/5)*5;

        % Draw the unity line 
        line([axis_min - 10, axis_max + 10],[axis_min - 10, axis_max + 10], ... 
            'Color', 'k', 'Linewidth', 1, 'Linestyle','--')

        for jj = 1:length(all_unit_names{xx,1})

            % If the unit significantly changed
            if sig_metric{xx,1}(jj,1) <= 0.05
                color_metric = sig_color;
            else % If the unit did not significantly change
                color_metric = insig_color;
            end

            % If the unit is a multi-unit
            if round(fract_contam{xx,1}(jj,1)) >= 0.1
                marker_metric = multi_marker;
                sz = multi_marker_size;
            else % If the unit is a single unit
                marker_metric = single_marker;
                sz = single_marker_size;
            end

            scatter(fire_rate_morn{xx}(jj), fire_rate_noon{xx}(jj), sz, marker_metric, 'MarkerEdgeColor', ... 
                [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);
            if isequal(unit_label, 1)
                text(fire_rate_morn{xx}(jj) + 1.5, fire_rate_noon{xx}(jj) - 1.5, ...
                    extractAfter(all_unit_names{xx}(jj), "elec"));
            end

            % Error
            err_morn = errorbar(fire_rate_morn{xx}(jj), fire_rate_noon{xx}(jj), ... 
                fire_rate_err_morn{xx}(jj), 'horizontal');
            err_morn.Color = [color_metric, 0, 0];
            err_morn.LineWidth = 1;
            err_noon = errorbar(fire_rate_morn{xx}(jj), fire_rate_noon{xx}(jj), ... 
                fire_rate_err_noon{xx}(jj), 'vertical');
            err_noon.Color = [color_metric, 0, 0];
            err_noon.LineWidth = 1;

        end % End of unit loop

        if strcmp(error_choice, 'best_fit')
            % Plot the standard deviation of the best fit line
            pos_best_fit_std = refline(units_slope_and_y_intercept(1), ... 
                (units_slope_and_y_intercept(2) + units_least_square_std));
            pos_best_fit_std.Color = 'r';
            pos_best_fit_std.LineStyle = '--';
            pos_best_fit_std.LineWidth = 1;
            neg_best_fit_std = refline(units_slope_and_y_intercept(1), ... 
                (units_slope_and_y_intercept(2) - units_least_square_std));
            neg_best_fit_std.Color = 'r';
            neg_best_fit_std.LineStyle = '--';
            neg_best_fit_std.LineWidth = 1;

            units_best_fit = refline(units_slope_and_y_intercept(1), units_slope_and_y_intercept(2));
            units_best_fit.Color = 'k';
            units_best_fit.LineWidth = 1;
        elseif strcmp(error_choice, 'ellip_err_prob') && length(fire_rate_morn{xx}) > 2
            % Plot the ellipital error
            for ii = 1:length(err_percent)
                plot(X_ellipse(ii,:),Y_ellipse(ii,:), 'Color', 'k', 'LineWidth',2);
            end
        end

        % Plot the hypothesis arrows
        if isequal(hypo_arrows, 1)
            ago_arrow_start = [axis_max - 10, axis_max];
            ago_arrow_stop = [axis_max, axis_max - 10];
            arrow(ago_arrow_start, ago_arrow_stop, 'Linewidth', 10, 'Length', 15, 'EdgeColor', [0 0.5 0]);
            antago_arrow_start = [axis_max - 5, axis_max - 15];
            antago_arrow_stop = [axis_max - 15, axis_max - 5];
            arrow(antago_arrow_start, antago_arrow_stop, 'Linewidth', 10, 'Length', 15, 'EdgeColor', 'r')
        end

        % Plot dummy points for the bottom right legend
        if isfield(Sampling_Params,'depth_sig') && strcmp(Sampling_Params.depth_sig, 'All')
            dummy_insig = plot(-50, -50, single_marker, 'MarkerSize', legend_size + 15, ...
                'MarkerEdgeColor', [insig_color, 0, 0], 'MarkerFaceColor',[1, 0, 0], 'LineWidth', 1.5);
            dummy_sig = plot(-55, -55, single_marker, 'MarkerSize', legend_size + 15, ...
                'MarkerEdgeColor', [sig_color, 0, 0], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);
        end

        % Plot dummy points for the top left legend
        dummy_single = plot(-57, -57, 'o', 'MarkerSize', legend_size - 15, ...
            'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5);
        if isfield(Sampling_Params,'ISI_quality') && strcmp(Sampling_Params.ISI_quality, 'All')
            dummy_multi = plot(-60, -60, multi_marker, 'MarkerSize', legend_size + 10, ...
                'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor', [1, 1, 1], 'LineWidth', 1.5);
        end

        % Set the axis
        xlim([axis_min - 10, axis_max + 10])
        ylim([axis_min - 10, axis_max + 10])

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
            title('');
            %title(fig_title, 'Fontsize', save_title_font_size - 10, 'Color', title_color);
        end

        % Plot the bottom right legend
        if strcmp(error_choice, 'best_fit')
            legend_r_square = strcat('r^2 =', {' '}, string(units_r_squared));
            legend_best_fit = units_best_fit;
            right_legend = legend([dummy_insig, dummy_sig, legend_best_fit], ...
                {'p > 0.05','p ≤ 0.05', legend_r_square}, ...
                'FontSize', legend_size, 'Location', 'southeast');
            right_legend.ItemTokenSize(1) = legend_size;
            legend boxoff
        else
            if isfield(Sampling_Params,'depth_sig') && strcmp(Sampling_Params.depth_sig, 'All')
                right_legend = legend([dummy_insig, dummy_sig], ...
                    {'p > 0.05','p ≤ 0.05'}, ...
                    'FontSize', legend_size, 'Location', 'southeast');
                right_legend.ItemTokenSize(1) = legend_size;
                legend boxoff
            end
        end
        
        % Only label every other tick
        x_labels = string(figure_axes.XAxis.TickLabels);
        y_labels = string(figure_axes.YAxis.TickLabels);
        x_labels(2:2:end) = NaN;
        y_labels(2:2:end) = NaN;
        figure_axes.XAxis.TickLabels = x_labels;
        figure_axes.YAxis.TickLabels = y_labels;

        % Reset the graph for the second legend
        legend_axes = axes('position', get(gca,'position'), 'visible', 'off');

        % Plot the top left legend
        if isfield(Sampling_Params,'ISI_quality') && strcmp(Sampling_Params.ISI_quality, 'All')
            left_legend = legend(legend_axes, [dummy_single, dummy_multi], {'Single Unit', 'Multi-Unit'}, ...
                'FontSize', legend_size, 'Location', 'northwest');
        else
            left_legend = legend(legend_axes, (dummy_single), {'Single Unit'}, ...
                'FontSize', legend_size, 'Location', 'northwest');
        end
        left_legend.ItemTokenSize(1) = legend_size;
        legend boxoff

        % Display the percent change & statistics
        fire_rate_mean_morn = mean(fire_rate_morn{xx,1});
        fire_rate_mean_noon = mean(fire_rate_noon{xx,1});
        [~, fire_rate_p_val] = ttest(fire_rate_morn{xx,1}, fire_rate_noon{xx,1});
        fire_rate_perc_change = ((fire_rate_mean_noon - fire_rate_mean_morn) / fire_rate_mean_morn) * 100;
        fprintf('Percent change is %0.2f%%, p = %0.3f \n', fire_rate_perc_change, fire_rate_p_val)

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







