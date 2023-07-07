function Anova_Plot(Monkey, Sampling_Params, Save_Figs)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('Baseline', 'Ramp', 'TgtHold', 'Peak', 'Depth')?
fire_rate_phase = 'Depth';

% Do you want to plot the morning / afternoon legend? (1 = Yes, 0 = No)
plot_legend = 1;

% What statistical measure do you want to use ('T-Test', 'Wilcox')
stat_test = 'T-Test';

% What effect size meausure do you want to use ('Perc_Change', 'Cohen')
effect_sz_test = 'Perc_Change';

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_Figs, 0)
    % Do you want a save title or blank title (1 = save_title, 0 = blank)
    Fig_Save_Title = 0;
    close all
end

%% Reassign variables according to what you're plotting

if strcmp(fire_rate_phase, 'Baseline')
    disp('Baseline Firing Rate')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'bsfr_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'bsfr_noon')};
end
if strcmp(fire_rate_phase, 'Ramp')
    disp('Ramp Phase')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'ramp_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'ramp_noon')};
end
if strcmp(fire_rate_phase, 'TgtHold')
    disp('TgtHold Phase')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'TgtHold_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'TgtHold_noon')};
end
if strcmp(fire_rate_phase, 'Peak')
    disp('Peak Firing Rate')
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
end
if strcmp(fire_rate_phase, 'Depth')
    disp('Depth of Modulation')
    fire_rate_morn = split_depth_excel{strcmp(column_names, 'depth_morn')};
    fire_rate_noon = split_depth_excel{strcmp(column_names, 'depth_noon')};
end

% Extract the other variables
drug_dose = split_depth_excel{strcmp(column_names, 'drug_dose_mg_per_kg')};
all_unit_names = split_depth_excel{strcmp(column_names, 'unit_names')};

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
fig_size = 700;

% Save Counter
ss = 1;

if ~isequal(Save_Figs, 0)
    save_title = strings;
end

%% Run through each experiment seperately

for xx = 1:length(all_unit_names)

    % Skip the function if no units match
    if isempty(all_unit_names{xx})
        disp('No units match criteria')
        continue
    end

    %% Add the monkey name to the title

    if strcmp(Monkey, 'All')
        fig_title = strcat('All Monkeys,', {' '});
    else
        fig_title = '';
        for ii = 1:length(Monkey)
            fig_title = strcat(fig_title, Monkey{ii}, ',', {' '});
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
    fire_rate_fig.Position = [200 50 fig_size fig_size];
    hold on

    histogram(fire_rate_morn{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
    histogram(fire_rate_noon{xx,1}, 15, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])

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
        'LineStyle','--', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', mean_line_width)
    line([fire_rate_mean_noon fire_rate_mean_noon], [y_limits(1) y_limits(2) + 0.25], ... 
        'LineStyle','--', 'Color', [.5 0 .5], 'LineWidth', mean_line_width)
        
    % Set the title
    if strcmp(Sampling_Params.trial_sessions, 'Ind')
        if ~strcmp(Sampling_Params.drug_choice, 'Con')
            title(strcat(fig_title, file_names{xx}, {' '}, string(drug_dose{xx}(1))), 'FontSize', title_font_size)
        else
            title(strcat(fig_title, file_names{xx}), 'FontSize', title_font_size)
        end
    elseif strcmp(Sampling_Params.trial_sessions, 'All')
        scatter_title = strcat('All Trials,', {' '}, Drug_Choice);
        if strcmp(Sampling_Params.trial_task, 'PG')
            scatter_title = strcat(scatter_title, {' '}, 'PG');
        end
        if strcmp(Sampling_Params.trial_task, 'KG')
            scatter_title = strcat(scatter_title, {' '}, 'KG');
        end
        if strcmp(Sampling_Params.trial_task, 'WS')
            scatter_title = strcat(scatter_title, {' '}, 'WS');
        end
        title(strcat(fig_title, scatter_title, {' '}, fire_rate_phase), 'FontSize', title_font_size)
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
    if strcmp(fire_rate_phase, 'Baseline')
        xlabel('Baseline Firing Rate (Hz)', 'FontSize', label_font_size)
    end
    if strcmp(fire_rate_phase, 'Ramp')
        xlabel('Ramp Phase Firing Rate (Hz)', 'FontSize', label_font_size)
    end
    if strcmp(fire_rate_phase, 'TgtHold')
        xlabel('TgtHold Phase Firing Rate (Hz)', 'FontSize', label_font_size)
    end
    if strcmp(fire_rate_phase, 'Peak')
        xlabel('Peak Firing Rate (Hz)', 'FontSize', label_font_size);
    end
    if strcmp(fire_rate_phase, 'Depth')
        xlabel('Depth of Modulation (Hz)', 'FontSize', label_font_size)
    end

    if ~isequal(Save_Figs, 0)
        fig_info = get(gca,'title');
        save_title(ss) = get(fig_info, 'string');
        % Make the title the drug
        if strcmp(Sampling_Params.drug_choice, 'Lex')
            title_color = [0 0.5 0];
            fig_title = strcat(fig_title, 'Escitalopram');
        end
        if strcmp(Sampling_Params.drug_choice, 'Caff')
            title_color = [0 0.5 0];
            fig_title = strcat(fig_title, 'Caffeine');
        end
        if strcmp(Sampling_Params.drug_choice, 'Cyp')
            title_color =  'r';
            fig_title = strcat(fig_title, 'Cyproheptadine');
        end
        if strcmp(Sampling_Params.drug_choice, 'Con')
            title_color =  'k';
            fig_title = strcat(fig_title, 'Control');
        end
        if contains(Sampling_Params.drug_choice, '202')
            title_color =  'k';
            fig_title = strcat(fig_title, file_names{xx});
        end
        if ~isequal(Fig_Save_Title, 0)
            if strcmp(Sampling_Params.trial_task, 'All')
                title(fig_title, 'FontSize', save_title_font_size - 15)
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


