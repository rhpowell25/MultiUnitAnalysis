function Plot_Unit_Metric(Monkey, Sampling_Params, Save_Figs)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which unit metric do you want to plot? 
% 'post_spike_facil', 'wave_sigtonoise', 'wave_peaktopeak', 
% 'alignment_morn', 'alignment_noon', 'fract_contam',
% 'nonlin_sigtonoise', 'nonlin_peaktopeak', 'spike_width', 'repol_time'
metric_choice = 'alignment_morn';

% Font specifications
label_font_size = 30;
legend_size = 25;
mean_line_width = 5;
n_value_dims = [0.55 0.45 0.44 0.44];
title_font_size = 14;
font_name = 'Arial';
fig_size = 700;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_Figs, 0)
    close all
end

%% Reassign variables according to what you're plotting

unit_metric = split_depth_excel{strcmp(column_names, metric_choice)};

% Extract the other variables
drug_dose = split_depth_excel{strcmp(column_names, 'drug_dose_mg_per_kg')};
all_unit_names = split_depth_excel{strcmp(column_names, 'unit_names')};

%% Run through each experiment seperately

for xx = 1:length(all_unit_names)

    % Skip the function if no units match
    if isempty(all_unit_names{xx})
        disp('No units match criteria')
        continue
    end

    fire_rate_n_val = length(all_unit_names{xx,1});

    %% Add the monkey name to the title

    if strcmp(Monkey, 'All')
        fig_title = strcat('All Monkeys,', {' '});
    else
        fig_title = '';
        for ii = 1:length(Monkey)
            fig_title = strcat(fig_title, Monkey{ii}, ',', {' '});
        end
    end

    %% Plot the unit metric histogram

    % Find the mean firing rates
    unit_metric_mean = mean(unit_metric{xx,1});

    % Merged
    unit_metric_fig = figure;
    unit_metric_fig.Position = [200 50 fig_size fig_size];
    hold on

    histogram(unit_metric{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [.7 .7 .7])

    % Set the axis
    x_limits = xlim;
    y_limits = ylim;
    xlim([x_limits(1), x_limits(2)])
    ylim([y_limits(1), y_limits(2) + 0.25])

    % Plot the mean
    line([unit_metric_mean, unit_metric_mean], [y_limits(1), y_limits(2) + 0.25], ... 
        'LineStyle','--', 'Color', 'k', 'LineWidth', mean_line_width)

    % Set the title
    if strcmp(Sampling_Params.trial_sessions, 'Ind')
        if ~strcmp(Sampling_Params.drug_choice, 'Con')
            title(strcat(fig_title, file_names{xx}, {' '}, string(drug_dose{xx}(1))), 'FontSize', title_font_size)
        else
            title(strcat(fig_title, file_names{xx}), 'FontSize', title_font_size)
        end
    elseif strcmp(Sampling_Params.trial_sessions, 'All')
        scatter_title = strcat('All Trials,', {' '}, Sampling_Params.drug_choice);
        if strcmp(Sampling_Params.trial_task, 'PG')
            scatter_title = strcat(scatter_title, {' '}, 'PG');
        end
        if strcmp(Sampling_Params.trial_task, 'KG')
            scatter_title = strcat(scatter_title, {' '}, 'KG');
        end
        if strcmp(Sampling_Params.trial_task, 'WS')
            scatter_title = strcat(scatter_title, {' '}, 'WS');
        end
        title(strcat(fig_title, scatter_title, {' '}, metric_choice), 'FontSize', title_font_size, 'Interpreter', 'none')
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
    xlabel(metric_choice, 'FontSize', label_font_size, 'Interpreter', 'none')
    
    % Annotation of the n_value
    n_value_string = strcat('n =', {' '}, mat2str(round(fire_rate_n_val, 3)));
    n_value_string = {char(n_value_string)};
    ann_n_value = annotation('textbox', n_value_dims, 'String', n_value_string, ...
        'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
        'EdgeColor','none', 'horizontalalignment', 'center');
    ann_n_value.FontSize = legend_size;
    ann_n_value.FontName = font_name;
    
    % Only label every other tick
    x_labels = string(figure_axes.XAxis.TickLabels);
    y_labels = string(figure_axes.YAxis.TickLabels);
    x_labels(2:2:end) = NaN;
    y_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    figure_axes.YAxis.TickLabels = y_labels;

end

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Documents\Grad School\Pop\20210922\XDS\Unsorted\PG\Nonlinear Energy Figures\';
    for ii = 1:numel(findobj('type','figure'))
        fig_info = get(gca,'title');
        fig_title = get(fig_info, 'string');
        fig_title = strrep(fig_title, ':', '');
        fig_title = strrep(fig_title, 'vs.', 'vs');
        fig_title = strrep(fig_title, 'mg.', 'mg');
        fig_title = strrep(fig_title, 'kg.', 'kg');
        fig_title = strrep(fig_title, '.', '_');
        fig_title = strrep(fig_title, '/', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(fig_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'png')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'pdf')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'fig')
        end
        close gcf
    end
end


