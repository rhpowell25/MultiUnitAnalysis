function Fire_Rate_Change_Histogram(Monkey, Sampling_Params, Save_File)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('Baseline', 'Ramp', 'TgtHold', 'Depth')?
fire_rate_phase = 'Depth';

% What statistical measure do you want to use ('T-Test', 'Wilcox')
stat_test = 'T-Test';

% How wide do you want the histogram edges to be
edge_width = 2;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_File, 0)
    close all
end

%% Reassign variables according to what you're plotting

if strcmp(fire_rate_phase, 'Baseline')
    disp('Baseline Firing Rate')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'bsfr_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'bsfr_noon')};
    fire_rate_t_test = split_depth_excel{strcmp(column_names, 'bsfr_t_test')};
    fire_rate_wilcoxon = split_depth_excel{strcmp(column_names, 'bsfr_wilcoxon')};
end
if strcmp(fire_rate_phase, 'Ramp')
    disp('Ramp Phase')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'ramp_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'ramp_noon')};
    fire_rate_t_test = split_depth_excel{strcmp(column_names, 'ramp_t_test')};
    fire_rate_wilcoxon = split_depth_excel{strcmp(column_names, 'ramp_wilcoxon')};
end
if strcmp(fire_rate_phase, 'TgtHold')
    disp('TgtHold Phase')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'TgtHold_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'TgtHold_noon')};
    fire_rate_t_test = split_depth_excel{strcmp(column_names, 'TgtHold_t_test')};
    fire_rate_wilcoxon = split_depth_excel{strcmp(column_names, 'TgtHold_wilcoxon')};
end
if strcmp(fire_rate_phase, 'Depth')
    disp('Depth of Modulation')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'depth_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'depth_noon')};
    fire_rate_t_test = split_depth_excel{strcmp(column_names, 'depth_t_test')};
    fire_rate_wilcoxon = split_depth_excel{strcmp(column_names, 'depth_wilcoxon')};
end

% Extract the other variables
drug_dose = split_depth_excel{strcmp(column_names, 'drug_dose_mg_per_kg')};
all_unit_names = split_depth_excel{strcmp(column_names, 'unit_names')};

%% Some variable extraction & definitions

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

% Bar color
if strcmp(Sampling_Params.drug_choice, 'Caff') || strcmp(Sampling_Params.drug_choice, 'Lex')
    bar_color = [0 0.5 0];
elseif strcmp(Sampling_Params.drug_choice, 'Cyp')
    bar_color = [1 0 0];
else
    bar_color = [0 0 0];
end

%% Loop through each of the experimental sessions
for xx = 1:length(all_unit_names)

    %% Add the monkey name to the title

    if strcmp(Monkey, 'All')
        Fig_Title = strcat('All Monkeys,', {' '});
    else
        Fig_Title = '';
        for ii = 1:length(Monkey)
            Fig_Title = strcat(Fig_Title, Monkey{ii}, ',', {' '});
        end
    end

    %% Do the statistics

    if strcmp(stat_test, 'T-Test')
        depth_p_value = fire_rate_t_test;
    elseif strcmp(stat_test, 'Wilcox')
        depth_p_value = fire_rate_wilcoxon;
    end

    %% Plot the multi_unit statistics graph

    fire_rate_change = fire_rate_noon{xx} - fire_rate_morn{xx};

    unit_stat_fig = figure;
    unit_stat_fig.Position = [200 50 800 600];
    hold on

    % Find the indices of the significant and insignificant depth changes
    sig_idx = find(depth_p_value{xx} <= 0.05);
    insig_idx = find(depth_p_value{xx} > 0.05);

    % Label the axis
    xlabel('Change In Firing Rate', 'FontSize', Plot_Params.label_font_size);
    ylabel('Units', 'FontSize', Plot_Params.label_font_size);

    % Set the edge size
    x_min = min(fire_rate_change);
    x_max = max(fire_rate_change);

    edge_min = round(x_min - 1);
    edge_max = round(x_max + 1);
    hist_edges = (edge_min:edge_width:edge_max);

    histogram(fire_rate_change(sig_idx), hist_edges, 'EdgeColor', 'k', 'FaceColor', bar_color)
    histogram(fire_rate_change(insig_idx), hist_edges, 'EdgeColor', 'k', 'FaceColor', [.7 .7 .7])

    % Calculate the mean change
    depth_mean = mean(fire_rate_change);
    sig_depth_mean = mean(fire_rate_change(sig_idx));
    insig_depth_mean = mean(fire_rate_change(insig_idx));
    fprintf("The mean change in modulation is %0.2f Hz \n", depth_mean);
    fprintf("The mean change of significant units is %0.2f Hz \n", sig_depth_mean);
    fprintf("The mean change of insignificant units is %0.2f Hz \n", insig_depth_mean);

    % Calculate the skew of the data
    depth_skew = skewness(fire_rate_change);
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
    figure_axes.FontSize = Plot_Params.label_font_size - 5;

    % Set the title
    if strcmp(Sampling_Params.trial_sessions, 'Ind')
        if ~strcmp(Sampling_Params.drug_choice, 'Con')
            title(strcat(Fig_Title, file_names{xx}, {' '}, string(drug_dose{xx}(1))), ...
                'FontSize', Plot_Params.title_font_size)
        else
            title(strcat(Fig_Title, file_names{xx}), 'FontSize', Plot_Params.title_font_size)
        end
    elseif strcmp(Sampling_Params.trial_sessions, 'All')
        scatter_title = strcat('All Tasks,', {' '}, Drug_Choice);
        if strcmp(Sampling_Params.trial_task, 'PG')
            scatter_title = strcat(scatter_title, {' '}, 'PG');
        end
        if strcmp(Sampling_Params.trial_task, 'WS')
            scatter_title = strcat(scatter_title, {' '}, 'WS');
        end
        title(strcat(Fig_Title, scatter_title, {' '}, 'Depth'), 'FontSize', Plot_Params.title_font_size)
    end

    % Draw the zero depth change line
    line([0,0],[y_limits(1), y_limits(2) + 1], 'Color', [.7 .7 .7], ...
        'Linewidth', Plot_Params.mean_line_width, 'Linestyle','--');
    % Draw the mean change in depth line
    line([depth_mean, depth_mean], [y_limits(1), y_limits(2) + 1], ... 
        'Color', bar_color, 'Linewidth', Plot_Params.mean_line_width - 1, 'Linestyle','--');

    % Plot dummy points for the right legend
    dummy_sig = plot(-12,-12, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', bar_color, 'MarkerFaceColor', bar_color, 'LineWidth', 1.5);
    dummy_grey = plot(-20,-20, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', [.7 .7 .7], 'MarkerFaceColor', [.7 .7 .7], 'LineWidth', 1.5);
 
    % Plot the right legend
    right_legend = legend([dummy_sig, dummy_grey], {'p ≤ 0.05', 'p > 0.05'}, ...
        'FontSize', Plot_Params.legend_size, 'Location', 'northeast');
    right_legend.ItemTokenSize(1) = Plot_Params.legend_size;
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

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end




