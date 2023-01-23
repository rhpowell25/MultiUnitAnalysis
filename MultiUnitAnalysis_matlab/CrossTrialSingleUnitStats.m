%% Load the output structures
% What Drug Do You Want To Load? ('Caff', 'Cyp', 'Lex', 'Con') or ('YYYYMMDD')
Drug_Choice = 'Con';
% What Monkey Do You Want To Load? ('All', 'Mihili', 'Jango', 'Jaco', 'Pop', 'Pancake')
Monkey_Choice = 'Pop';
% What targets do you want to use ('Max', 'Min')
Tgt_Choice = 'Max';

[~, ~, depth_morn, depth_noon, ~, ~, ~, ~, ~, depth_p_value, ~, nonlin_p_value, fract_contam, ...
    pref_dir, ~, drug_dose, file_name, all_unit_names] = Fast_Load_Depth(Drug_Choice, Monkey_Choice, Tgt_Choice);

%% Display the function being used
disp('Cross Trial Single Unit Statistics:')

%% Some of the plotting specifications
% Individual Sessions or All Sessions? ('per_trial' vs 'all_trials')
trial_sessions = 'per_trial';

% Unit quality ('All' vs. 'Stable')
unit_quality = 'Stable';

% ISI quality ('All' vs. 'Single')
ISI_quality = 'All';

% Which task do you want to plot ('PG' vs. 'WS', vs 'Both)?
trial_task = 'Both';

% What preferred direction do you want to plot(-90, 0, 90, 180, 'All')
plot_dir = 'All';

% What minimum depth of modulation do you want to observe
depth_min = NaN;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
Save_Figs = 'pdf';
if ~isequal(Save_Figs, 0)
    % Do you want a save title or blank title (1 = save_title, 0 = blank)
    Fig_Save_Title = 0;
    close all
end

% How wide do you want the histogram edges to be
edge_width = 2;

%% Some variable extraction & definitions

% Font specifications
label_font_size = 30;
legend_size = 20;
mean_line_width = 5;
title_font_size = 14;
save_title_font_size = 45;

% Save Counter
ss = 1;

% Bar color
if strcmp(Drug_Choice, 'Caff') || strcmp(Drug_Choice, 'Lex')
    bar_color = [0 0.5 0];
elseif strcmp(Drug_Choice, 'Cyp')
    bar_color = [1 0 0];
else
    bar_color = [0 0 0];
end

if ~isequal(Save_Figs, 0)
    save_title = strings;
end

%% Only use trial tasks selected
if strcmp(trial_task, 'PG')
    disp('Powergrasp Only')
    powergrasp_idx = find(~contains(file_name, 'PG'));
    all_unit_names(powergrasp_idx) = [];
    depth_morn(powergrasp_idx) = [];
    depth_noon(powergrasp_idx) = [];
    depth_p_value(powergrasp_idx) = [];
    nonlin_p_value(powergrasp_idx) = [];
    fract_contam(powergrasp_idx) = [];
    pref_dir(powergrasp_idx) = [];
    file_name(powergrasp_idx) = [];
end

if strcmp(trial_task, 'WS')
    disp('Wrist w/ Spring Only')
    wrist_idx = find(~contains(file_name, 'WS'));
    all_unit_names(wrist_idx) = [];
    depth_morn(wrist_idx) = [];
    depth_noon(wrist_idx) = [];
    depth_p_value(wrist_idx) = [];
    nonlin_p_value(wrist_idx) = [];
    fract_contam(wrist_idx) = [];
    pref_dir(wrist_idx) = [];
    file_name(wrist_idx) = [];
end

%% Only use the preferred direction selected

if ~strcmp(plot_dir, 'All')
    fprintf('Preferred Direciton of %.0f° Only \n', plot_dir)
    for ii = 1:length(pref_dir)
        pref_dir_idx = find(pref_dir{ii,1} ~= plot_dir);
        all_unit_names{ii,1}(pref_dir_idx) = [];
        depth_morn{ii,1}(pref_dir_idx) = [];
        depth_noon{ii,1}(pref_dir_idx) = [];
        depth_p_value{ii,1}(pref_dir_idx) = [];
        nonlin_p_value{ii,1}(pref_dir_idx) = [];
        fract_contam{ii,1}(pref_dir_idx) = [];
    end
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
    all_trials_depth_morn = zeros(length(total_units),1);
    all_trials_depth_noon = zeros(length(total_units),1);
    all_trials_depth_p_value = zeros(length(total_units),1);
    all_trials_sort_p_value = zeros(length(total_units),1);
    all_trials_merged_ISI_mode = zeros(length(total_units),1);
    cc = 1;
    for xx = 1:length(all_unit_names)
        for jj = 1:length(all_unit_names{xx})
            all_trials_unit_names{cc,1} = all_unit_names{xx}(jj);
            all_trials_depth_morn(cc,1) = depth_morn{xx}(jj);
            all_trials_depth_noon(cc,1) = depth_noon{xx}(jj);
            all_trials_depth_p_value(cc,1) = depth_p_value{xx}(jj);
            all_trials_sort_p_value(cc,1) = nonlin_p_value{xx}(jj);
            all_trials_merged_ISI_mode(cc,1) = fract_contam{xx}(jj);
            cc = cc + 1;
        end
    end
    % Rename the now merged variables
    all_unit_names = struct([]);
    all_unit_names{1,1} = all_trials_unit_names;
    depth_morn = struct([]);
    depth_morn{1,1} = all_trials_depth_morn;
    depth_noon = struct([]);
    depth_noon{1,1} = all_trials_depth_noon;
    depth_p_value = struct([]);
    depth_p_value{1,1} = all_trials_depth_p_value;
    nonlin_p_value = struct([]);
    nonlin_p_value{1,1} = all_trials_sort_p_value;
    fract_contam = struct([]);
    fract_contam{1,1} = all_trials_merged_ISI_mode;
end

%% Remove any units that fall below the minimum selected depth of modulation

if ~isnan(depth_min)
    min_depth_violations_morn = struct([]);
    min_depth_violations_noon = struct([]);
    min_depth_violations = struct([]);
    for ii = 1:length(all_unit_names)
        min_depth_violations_morn{ii,1} = find(depth_morn{ii,1} < depth_min);
        min_depth_violations_noon{ii,1} = find(depth_noon{ii,1} < depth_min);
        min_depth_violations{ii,1} = intersect(min_depth_violations_morn{ii,1}, min_depth_violations_noon{ii,1});
        all_unit_names{ii,1}(min_depth_violations{ii,1}) = [];
        depth_morn{ii,1}(min_depth_violations{ii,1}) = [];
        depth_noon{ii,1}(min_depth_violations{ii,1}) = [];
        nonlin_p_value{ii,1}(min_depth_violations{ii,1}) = [];
        fract_contam{ii,1}(min_depth_violations{ii,1}) = [];
    end
end

%% Loop through each of the experimental sessions
for xx = 1:length(all_unit_names)

    %% Find the indexes of the poorly sorted units
    poor_fract_contam_idx = find(fract_contam{xx,1} >= 0.1);
    sig_p_value_idx = find(nonlin_p_value{xx} <= 0.05);
    all_units_idx = 1:length(all_unit_names{xx,1});
    poorly_sorted_idx = unique(cat(1, poor_fract_contam_idx, sig_p_value_idx));

    all_single_units_idx = all_units_idx;
    all_single_units_idx(poor_fract_contam_idx) = [];

    all_best_units_idx = all_units_idx;
    all_best_units_idx(sig_p_value_idx) = [];

    best_single_units_idx = all_units_idx;
    best_single_units_idx(poorly_sorted_idx) = [];
    
    if strcmp(unit_quality, 'All') && strcmp(ISI_quality, 'All')
        units_idx = all_units_idx;
        units_title = strcat('All Units', {' '});
    end

    if strcmp(unit_quality, 'All') && strcmp(ISI_quality, 'Single')
        units_idx = all_single_units_idx;
        units_title = strcat('Single Units,', {' '});
    end

    if strcmp(unit_quality, 'Stable') && strcmp(ISI_quality, 'All')
        units_idx = all_best_units_idx;
        units_title = strcat('Stable Units,', {' '});
    end

    if strcmp(unit_quality, 'Stable') && strcmp(ISI_quality, 'Single')
        units_idx = best_single_units_idx;
        units_title = strcat('Stable, Single Units,', {' '});
    end

    %% Plot the multi_unit statistics graph

    depth_change = depth_noon{xx} - depth_morn{xx};

    unit_stat_fig = figure;
    unit_stat_fig.Position = [200 50 800 600];
    hold on

    % Set the title
    if strcmp(trial_sessions, 'per_trial')
        title(strcat(units_title, file_name{xx}, {' '}, drug_dose(xx)), 'FontSize', title_font_size)
    elseif strcmp(trial_sessions, 'all_trials')
        scatter_title = 'All Trials';
        title(strcat(units_title, scatter_title, {' '}, Drug_Choice), 'FontSize', title_font_size)
    end

    % Find the indices of the significant and insignificant depth changes
    sig_idx = find(depth_p_value{xx}(units_idx) <= 0.05);
    insig_idx = find(depth_p_value{xx}(units_idx) > 0.05);

    % Label the axis
    xlabel('Change In Depth of Modulation', 'FontSize', label_font_size);
    ylabel('Units', 'FontSize', label_font_size);

    % Set the edge size
    x_min = min(depth_change);
    x_max = max(depth_change);

    edge_min = round(x_min - 1);
    edge_max = round(x_max + 1);
    hist_edges = (edge_min:edge_width:edge_max);

    histogram(depth_change(units_idx(sig_idx)), hist_edges, 'EdgeColor', 'k', 'FaceColor', bar_color)
    histogram(depth_change(units_idx(insig_idx)), hist_edges, 'EdgeColor', 'k', 'FaceColor', [.7 .7 .7])

    % Calculate the mean change
    depth_mean = mean(depth_change(units_idx));
    sig_depth_mean = mean(depth_change(units_idx(sig_idx)));
    insig_depth_mean = mean(depth_change(units_idx(insig_idx)));
    fprintf("The mean change in modulation is %0.2f Hz \n", depth_mean);
    fprintf("The mean change of significant units is %0.2f Hz \n", sig_depth_mean);
    fprintf("The mean change of insignificant units is %0.2f Hz \n", insig_depth_mean);

    % Calculate the skew of the data
    depth_skew = skewness(depth_change(units_idx));
    fprintf("The skewness of the distribution is %0.2f \n", depth_skew);

    % Set the axis
    x_limits = max(abs(x_min), abs(x_max)) + 5;
    xlim([-x_limits, x_limits])
    y_limits = ylim;
    ylim([y_limits(1), y_limits(2)])

    % Only label every other tick
    figure_axes = gca;
    % Set ticks to outside
    set(figure_axes,'TickDir','out');
    % Remove the top and right tick marks
    set(figure_axes,'box','off')
    % Set the tick label font size
    figure_axes.FontSize = label_font_size - 5;

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
            if strcmp(trial_task, 'Both')
                title(strcat(units_title, fig_title), 'FontSize', save_title_font_size - 15)
            else
                if strcmp(plot_dir, 'Both')
                    title(strcat(units_title, fig_title, {' '}, trial_task), 'FontSize', save_title_font_size - 15)
                else
                    title(strcat(units_title, fig_title, {' '}, trial_task, {' '}, num2str(plot_dir), '°'), ...
                        'FontSize', save_title_font_size - 15)
                    save_title(ss) = get(fig_info, 'string');
                end
            end
            title(fig_title, 'Fontsize', save_title_font_size, 'Color', bar_color);
        else
            title('')
        end
    end

    % Draw the zero depth change line
    zero_depth = line([0,0],[y_limits(1), y_limits(2) + 1], 'Color', [.7 .7 .7], 'Linewidth', mean_line_width, 'Linestyle','--');
    % Draw the mean change in depth line
    mean_depth = line([depth_mean, depth_mean], [y_limits(1), y_limits(2) + 1], ... 
        'Color', bar_color, 'Linewidth', mean_line_width - 1, 'Linestyle','--');

    % Plot dummy points for the right legend
    dummy_sig = plot(-12,-12, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', bar_color, 'MarkerFaceColor', bar_color, 'LineWidth', 1.5);
    dummy_grey = plot(-20,-20, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', [.7 .7 .7], 'MarkerFaceColor', [.7 .7 .7], 'LineWidth', 1.5);
 
    % Plot the right legend
    right_legend = legend([dummy_sig, dummy_grey], {'p ≤ 0.05', 'p > 0.05'}, ...
        'FontSize', legend_size, 'Location', 'northeast');
    right_legend.ItemTokenSize(1) = legend_size;
    legend boxoff

    % Convert the Y-ticks to percent
    ytix = get(figure_axes, 'YTick');
    set(gca, 'YTick',ytix, 'YTickLabel',ytix)
    x_labels = string(figure_axes.XAxis.TickLabels);
    y_labels = string(figure_axes.YAxis.TickLabels);
    x_labels(2:2:end) = NaN;
    y_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    figure_axes.YAxis.TickLabels = y_labels;

    % Add to the counter
    ss = ss + 1;

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
        save_title(ii) = strcat(save_title(ii), {' '}, '(Stats)');
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



