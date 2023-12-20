function Plot_Unit_Metric(Monkey, Sampling_Params, Save_File)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which unit metric do you want to plot? 
% 'post_spike_facil', 'wave_sigtonoise', 'wave_peaktopeak', 
% 'alignment_morn', 'alignment_noon', 'fract_contam',
% 'nonlin_sigtonoise', 'nonlin_peaktopeak', 'spike_width', 'repol_time'
metric_choice = 'alignment_morn';

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;
n_value_dims = [0.55 0.45 0.44 0.44];

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_File, 0)
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
        Fig_Title = strcat('All Monkeys,', {' '});
    else
        Fig_Title = '';
        for ii = 1:length(Monkey)
            Fig_Title = strcat(Fig_Title, Monkey{ii}, ',', {' '});
        end
    end

    %% Plot the unit metric histogram

    % Find the mean firing rates
    unit_metric_mean = mean(unit_metric{xx,1});

    % Merged
    unit_metric_fig = figure;
    unit_metric_fig.Position = [200 50 Plot_Params.fig_size Plot_Params.fig_size];
    hold on

    histogram(unit_metric{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [.7 .7 .7])

    % Set the axis
    x_limits = xlim;
    y_limits = ylim;
    xlim([x_limits(1), x_limits(2)])
    ylim([y_limits(1), y_limits(2) + 0.25])

    % Plot the mean
    line([unit_metric_mean, unit_metric_mean], [y_limits(1), y_limits(2) + 0.25], ... 
        'LineStyle','--', 'Color', 'k', 'LineWidth', Plot_Params.mean_line_width)

    % Set the title
    if strcmp(Sampling_Params.trial_sessions, 'Ind')
        if ~strcmp(Sampling_Params.drug_choice, 'Con')
            title(strcat(Fig_Title, file_names{xx}, {' '}, string(drug_dose{xx}(1))), ...
                'FontSize', Plot_Params.title_font_size)
        else
            title(strcat(Fig_Title, file_names{xx}), 'FontSize', Plot_Params.title_font_size)
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
        title(strcat(Fig_Title, scatter_title, {' '}, metric_choice), ...
            'FontSize', Plot_Params.title_font_size, 'Interpreter', 'none')
    end

    % Axis Editing
    figure_axes = gca;
    % Set ticks to outside
    set(figure_axes,'TickDir','out');
    % Remove the top and right tick marks
    set(figure_axes,'box','off')
    % Set the tick label font size
    figure_axes.FontSize = Plot_Params.label_font_size - 15;

    % Label the axis
    ylabel('Units', 'FontSize', Plot_Params.label_font_size)
    xlabel(metric_choice, 'FontSize', Plot_Params.label_font_size, 'Interpreter', 'none')
    
    % Annotation of the n_value
    n_value_string = strcat('n =', {' '}, mat2str(round(fire_rate_n_val, 3)));
    n_value_string = {char(n_value_string)};
    ann_n_value = annotation('textbox', n_value_dims, 'String', n_value_string, ...
        'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
        'EdgeColor','none', 'horizontalalignment', 'center');
    ann_n_value.FontSize = Plot_Params.legend_size;
    ann_n_value.FontName = Plot_Params.font_name;
    
    % Only label every other tick
    x_labels = string(figure_axes.XAxis.TickLabels);
    y_labels = string(figure_axes.YAxis.TickLabels);
    x_labels(2:2:end) = NaN;
    y_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    figure_axes.YAxis.TickLabels = y_labels;

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end


