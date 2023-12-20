function Histogram_Plot(Monkey, Sampling_Params, Save_File)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('bsfr', 'ramp', 'TgtHold', 'Peak', 'depth')?
fire_rate_phase = 'depth';

% Do you want to plot the morning / afternoon legend? (1 = Yes, 0 = No)
plot_legend = 1;

% What statistical measure do you want to use ('T-Test', 'Wilcox')
stat_test = 'T-Test';

% What effect size meausure do you want to use ('Perc_Change', 'Cohen')
effect_sz_test = 'Perc_Change';

% Do you want the simplified title? (1 = Yes, 0 = No)
simp_title = 0;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_File, 0)
    close all
end

%% Reassign variables according to what you're plotting

if strcmp(fire_rate_phase, 'Peak')
    bsfr_morn = split_depth_excel{strcmp(column_names, 'bsfr_morn')};
    bsfr_noon = split_depth_excel{strcmp(column_names, 'bsfr_noon')};
    depth_morn = split_depth_excel{strcmp(column_names, 'depth_morn')};
    depth_noon = split_depth_excel{strcmp(column_names, 'depth_noon')};
    fire_rate_morn = struct([]);
    fire_rate_noon = struct([]);
    for ii = 1:length(bsfr_morn)
        fire_rate_morn{ii} = bsfr_morn{ii} + depth_morn{ii};
        fire_rate_noon{ii} = bsfr_noon{ii} + depth_noon{ii};
    end
else
    fire_rate_morn = split_depth_excel{strcmp(column_names, strcat(fire_rate_phase, '_morn'))};
    fire_rate_noon = split_depth_excel{strcmp(column_names, strcat(fire_rate_phase, '_noon'))};
end

% Extract the other variables
drug_dose = split_depth_excel{strcmp(column_names, 'drug_dose_mg_per_kg')};
all_unit_names = split_depth_excel{strcmp(column_names, 'unit_names')};

%% Some variable extraction & definitions

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;
if isequal(plot_legend, 1)
    p_value_dims = [0.51 0.3 0.44 0.44];
    n_value_dims = [0.49 0.225 0.44 0.44];
else
    p_value_dims = [0.51 0.45 0.44 0.44];
    n_value_dims = [0.49 0.375 0.44 0.44];
end

%% Run through each experiment seperately

for xx = 1:length(all_unit_names)

    % Skip the function if no units match
    if isempty(all_unit_names{xx})
        disp('No units match criteria')
        continue
    end

    %% Axis & Title information

    % Title info
    Fig_Title = '';
    if strcmp(Monkey, 'All')
        Fig_Title = strcat('All Monkeys,', {' '});
    end
    if strcmp(Sampling_Params.trial_sessions, 'All')
        for ii = 1:length(Monkey)
            Fig_Title = strcat(Fig_Title, Monkey{ii}, ',', {' '});
        end
        Fig_Title = strcat(Fig_Title, {' '}, 'All Trials,', {' '}, Sampling_Params.drug_choice);
        if strcmp(Sampling_Params.trial_task, 'PG')
            Fig_Title = strcat(Fig_Title, {' '}, 'PG');
        end
        if strcmp(Sampling_Params.trial_task, 'KG')
            Fig_Title = strcat(Fig_Title, {' '}, 'KG');
        end
        if strcmp(Sampling_Params.trial_task, 'WS')
            Fig_Title = strcat(Fig_Title, {' '}, 'WS');
        end
        Fig_Title = strcat(Fig_Title, {' '}, fire_rate_phase);
    else
        if ~strcmp(Sampling_Params.drug_choice, 'Con')
            Fig_Title = strcat(Fig_Title, file_names{xx}, {' '}, string(drug_dose{xx}(1)));
        else
            Fig_Title = strcat(Fig_Title, file_names{xx});
        end
    end
    
    % Simplified title
    if ~isequal(simp_title, 0)
        if strcmp(Sampling_Params.drug_choice, 'Lex')
            Plot_Params.title_color = [0 0.5 0];
            Fig_Title = 'Escitalopram';
        end
        if strcmp(Sampling_Params.drug_choice, 'Caff')
            Plot_Params.title_color = [0 0.5 0];
            Fig_Title = 'Caffeine';
        end
        if strcmp(Sampling_Params.drug_choice, 'Cyp')
            Plot_Params.title_color =  'r';
            Fig_Title = 'Cyproheptadine';
        end
        if strcmp(Sampling_Params.drug_choice, 'Con')
            Plot_Params.title_color =  'k';
            Fig_Title = 'Control';
        end
        if contains(Sampling_Params.drug_choice, '202')
            Fig_Title = strcat(Fig_Title, file_names{xx});
        end
    end

    % Label info
    if strcmp(fire_rate_phase, 'bsfr')
        x_label = 'Baseline Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'ramp')
        x_label = 'Ramp Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'TgtHold')
        x_label = 'TgtHold Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'Peak')
        x_label = 'Peak Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'depth')
        x_label = 'Depth of Modulation (Hz)';
    end
    if strcmp(fire_rate_phase, 'EMG_amp')
        x_label = 'Peak EMG Amplitude';
    end

    %% Do the statistics

    if strcmp(stat_test, 'T-Test')
        [~, fire_rate_p_val] = ttest(fire_rate_morn{xx,1}, fire_rate_noon{xx,1});
    elseif strcmp(stat_test, 'Wilcox')
        [fire_rate_p_val, ~] = ranksum(fire_rate_morn{xx,1}, fire_rate_noon{xx,1});
    end
    fire_rate_n_val = length(fire_rate_morn{xx,1});

    %% Plot the depth of modulation histograms

    % Find the mean firing rates
    fire_rate_mean_morn = mean(fire_rate_morn{xx,1});
    fire_rate_mean_noon = mean(fire_rate_noon{xx,1});
    if strcmp(effect_sz_test, 'Perc_Change')
        fire_rate_effect_size = ((fire_rate_mean_noon - fire_rate_mean_morn) / fire_rate_mean_morn) * 100;
    elseif strcmp(effect_sz_test, 'Cohen')
        fire_rate_effect_size = Cohen_D(fire_rate_morn{xx,1}, fire_rate_noon{xx,1});
    end

    % Merged
    fire_rate_fig = figure;
    fire_rate_fig.Position = [200 50 Plot_Params.fig_size Plot_Params.fig_size];
    hold on

    histogram(fire_rate_morn{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
    histogram(fire_rate_noon{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])

    disp(strcat(Monkey, {' '}, Sampling_Params.drug_choice, ':'))
    fprintf('Effect size is %0.2f, p = %0.3f \n', fire_rate_effect_size, fire_rate_p_val)

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
        'LineStyle','--', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', Plot_Params.mean_line_width)
    line([fire_rate_mean_noon fire_rate_mean_noon], [y_limits(1) y_limits(2) + 0.25], ... 
        'LineStyle','--', 'Color', [.5 0 .5], 'LineWidth', Plot_Params.mean_line_width)
        
    % Set the title
    title(Fig_Title, 'FontSize', Plot_Params.title_font_size, 'Color', ...
        Plot_Params.title_color, 'Interpreter', 'none')
    %title('')

    % Axis Editing
    figure_axes = gca;
    % Set ticks to outside
    set(figure_axes,'TickDir','out');
    % Remove the top and right tick marks
    set(figure_axes,'box','off')
    % Set the tick label font size
    figure_axes.FontSize = Plot_Params.label_font_size;
    % Set The Font
    set(figure_axes,'fontname', Plot_Params.font_name);

    % Label the axis
    ylabel('Units', 'FontSize', Plot_Params.label_font_size)
    xlabel(x_label, 'FontSize', Plot_Params.label_font_size);

    % Annotation of the p_value
    if round(fire_rate_p_val, 3) > 0
        p_value_string = strcat('p =', {' '}, mat2str(round(fire_rate_p_val, 3)));
        p_value_string = {char(p_value_string)};
        ann_p_value = annotation('textbox', p_value_dims, 'String', p_value_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_p_value.FontSize = Plot_Params.legend_size;
        ann_p_value.FontName = Plot_Params.font_name;
    end

    if isequal(round(fire_rate_p_val, 3), 0)
        p_value_string = strcat('p <', {' '}, '0.001');
        p_value_string = {char(p_value_string)};
        ann_p_value = annotation('textbox', p_value_dims, 'String', p_value_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_p_value.FontSize = Plot_Params.legend_size;
        ann_p_value.FontName = Plot_Params.font_name;
    end

    % Annotation of the n_value
    n_value_string = strcat('n =', {' '}, mat2str(round(fire_rate_n_val, 3)));
    n_value_string = {char(n_value_string)};
    ann_n_value = annotation('textbox', n_value_dims, 'String', n_value_string, ...
        'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
        'EdgeColor','none', 'horizontalalignment', 'center');
    ann_n_value.FontSize = Plot_Params.legend_size;
    ann_n_value.FontName = Plot_Params.font_name;

    if isequal(plot_legend, 1)
        legend([dummy_morn, dummy_noon], ... 
            {'Morning', 'Afternoon'}, ... 
            'FontSize', Plot_Params.legend_size, 'Location', 'NorthEast')
        legend boxoff
    end

    % Only label every other tick
    x_labels = string(figure_axes.XAxis.TickLabels);
    y_labels = string(figure_axes.YAxis.TickLabels);
    x_labels(2:2:end) = NaN;
    y_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    figure_axes.YAxis.TickLabels = y_labels;

    %% Save the file if selected
    Fig_Title = strcat(Fig_Title, {' '}, 'Hist');
    Save_Figs(Fig_Title, Save_File)

end

