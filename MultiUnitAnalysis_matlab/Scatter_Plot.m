function Scatter_Plot(Monkey, Sampling_Params, Save_File)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('bsfr', 'ramp', 'TgtHold', 'Peak', 'depth', 'EMG_amp')?
fire_rate_phase = 'depth';

% Do you want to plot the legends? (1 = Yes, 0 = No)
plot_legends = 1;

% Best fit or elliptial error probable? ('none', 'best_fit', or 'ellip_err_prob')
error_choice = 'none';

% What effect size meausure do you want to use ('Perc_Change', 'Cohen')
effect_sz_test = 'Perc_Change';

% Do you want the name of each unit labeled? (1 = Yes, 0 = No)
unit_label = 0;

% Do you want to plot the hypothesis arrows? (1 = Yes, 0 = No)
hypo_arrows = 0;

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
    fire_rate_err_morn = split_depth_excel{strcmp(column_names, 'depth_err_morn')};
    fire_rate_err_noon = split_depth_excel{strcmp(column_names, 'depth_err_noon')};
    fire_rate_p_val = split_depth_excel{strcmp(column_names, 'depth_p_val')};
else
    fire_rate_morn = split_depth_excel{strcmp(column_names, strcat(fire_rate_phase, '_morn'))};
    fire_rate_err_morn = split_depth_excel{strcmp(column_names, strcat(fire_rate_phase, '_err_morn'))};
    fire_rate_noon = split_depth_excel{strcmp(column_names, strcat(fire_rate_phase, '_noon'))};
    fire_rate_err_noon = split_depth_excel{strcmp(column_names, strcat(fire_rate_phase, '_err_noon'))};
    fire_rate_p_val = split_depth_excel{strcmp(column_names, strcat(fire_rate_phase, '_p_val'))};
end

% Extract the other variables
fract_contam = split_depth_excel{strcmp(column_names, 'fract_contam')};
drug_dose = split_depth_excel{strcmp(column_names, 'drug_dose_mg_per_kg')};
all_unit_names = split_depth_excel{strcmp(column_names, 'unit_names')};
if contains('EMG_names', column_names)
    all_emg_names = split_depth_excel{strcmp(column_names, 'EMG_names')};
end

%% Some variable extraction & definitions

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

% Scatter Marker Shapes
single_marker = '.';
multi_marker = '*';

% Scatter Marker Colors
insig_color = 0;
sig_color = 1;

% Scatter Marker sizes
single_marker_size = 1000;
multi_marker_size = 250;

%% Loop through each of the experimental sessions
for xx = 1:length(all_unit_names)

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
        label_addition = 'Baseline Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'ramp')
        label_addition = 'Ramp Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'TgtHold')
        label_addition = 'TgtHold Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'Peak')
        label_addition = 'Peak Firing Rate (Hz)';
    end
    if strcmp(fire_rate_phase, 'depth')
        label_addition = 'Depth of Modulation (Hz)';
    end
    if strcmp(fire_rate_phase, 'EMG_amp')
        label_addition = 'Peak EMG Amplitude';
    end

    %% Plot the Depth of Modulation Scatter

    scatter_fig = figure;
    scatter_fig.Position = [200 50 Plot_Params.fig_size Plot_Params.fig_size];
    hold on

    % Plot a circle (x-a) = rcos(t) / (y-b) = rsin(t)
    %t = linspace(0,2*pi);
    %rad = 7;
    %x = rad*cos(t) + fire_rate_morn{1,1}(12);
    %y = rad*sin(t) + fire_rate_noon{1,1}(12);
    %plot(x,y, 'r', 'LineWidth', 2)

    % Set the title
    title(Fig_Title, 'FontSize', Plot_Params.title_font_size, ...
        'Color', Plot_Params.title_color, 'Interpreter', 'none')
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
    xlabel(strcat('Morning', {' '}, label_addition), 'FontSize', Plot_Params.label_font_size);
    ylabel(strcat('Afternoon', {' '}, label_addition), 'FontSize', Plot_Params.label_font_size);

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
        if fire_rate_p_val{xx,1}(jj,1) <= 0.05
            color_metric = sig_color;
        else % If the unit did not significantly change
            color_metric = insig_color;
        end

        % If the unit is a multi-unit
        if fract_contam{xx,1}(jj,1) >= 0.1
            marker_metric = multi_marker;
            sz = multi_marker_size;
        else % If the unit is a single unit
            marker_metric = single_marker;
            sz = single_marker_size;
        end
        if strcmp(fire_rate_phase, 'EMG_amp')
            marker_metric = single_marker;
            sz = single_marker_size;
        end

        scatter(fire_rate_morn{xx}(jj), fire_rate_noon{xx}(jj), sz, marker_metric, 'MarkerEdgeColor', ... 
            [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);
        if isequal(unit_label, 1)
            if strcmp(fire_rate_phase, 'EMG_amp')
                text(fire_rate_morn{xx}(jj) + 1.5, fire_rate_noon{xx}(jj) - 1.5, ...
                    all_emg_names{xx}(jj));
            else
                text(fire_rate_morn{xx}(jj) + 1.5, fire_rate_noon{xx}(jj) - 1.5, ...
                    extractAfter(all_unit_names{xx}(jj), "elec"));
            end
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
    if isequal(plot_legends, 1)
        if isfield(Sampling_Params,'depth_sig') && strcmp(Sampling_Params.depth_sig, 'All')
            dummy_insig = plot(-50, -50, single_marker, 'MarkerSize', Plot_Params.legend_size + 15, ...
                'MarkerEdgeColor', [insig_color, 0, 0], 'MarkerFaceColor',[1, 0, 0], 'LineWidth', 1.5);
            dummy_sig = plot(-55, -55, single_marker, 'MarkerSize', Plot_Params.legend_size + 15, ...
                'MarkerEdgeColor', [sig_color, 0, 0], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);
        end
        % Plot dummy points for the top left legend
        dummy_single = plot(-57, -57, 'o', 'MarkerSize', Plot_Params.legend_size - 15, ...
            'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'LineWidth', 1.5);
        if isfield(Sampling_Params,'ISI_quality') && strcmp(Sampling_Params.ISI_quality, 'All')
            dummy_multi = plot(-60, -60, multi_marker, 'MarkerSize', Plot_Params.legend_size + 10, ...
                'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor', [1, 1, 1], 'LineWidth', 1.5);
        end
    end

    % Set the axis
    xlim([axis_min - 10, axis_max + 10])
    ylim([axis_min - 10, axis_max + 10])

    % Plot the bottom right legend
    if isequal(plot_legends, 1)
        if strcmp(error_choice, 'best_fit')
            legend_r_square = strcat('r^2 =', {' '}, string(units_r_squared));
            legend_best_fit = units_best_fit;
            right_legend = legend([dummy_insig, dummy_sig, legend_best_fit], ...
                {'p > 0.05','p ≤ 0.05', legend_r_square}, ...
                'FontSize', Plot_Params.legend_size, 'Location', 'southeast');
            right_legend.ItemTokenSize(1) = Plot_Params.legend_size;
            legend boxoff
        else
            if isfield(Sampling_Params,'depth_sig') && strcmp(Sampling_Params.depth_sig, 'All')
                right_legend = legend([dummy_insig, dummy_sig], ...
                    {'p > 0.05','p ≤ 0.05'}, ...
                    'FontSize', Plot_Params.legend_size, 'Location', 'southeast');
                right_legend.ItemTokenSize(1) = Plot_Params.legend_size;
                legend boxoff
            end
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
    if isequal(plot_legends, 1)
        legend_axes = axes('position', get(gca,'position'), 'visible', 'off');
    end

    % Plot the top left legend
    if isequal(plot_legends, 1)
        if isfield(Sampling_Params,'ISI_quality') && strcmp(Sampling_Params.ISI_quality, 'All')
            left_legend = legend(legend_axes, [dummy_single, dummy_multi], {'Single Unit', 'Multi-Unit'}, ...
                'FontSize', Plot_Params.legend_size, 'Location', 'northwest');
        else
            left_legend = legend(legend_axes, (dummy_single), {'Single Unit'}, ...
                'FontSize', Plot_Params.legend_size, 'Location', 'northwest');
        end
        left_legend.ItemTokenSize(1) = Plot_Params.legend_size;
        legend boxoff
    end

    % Display the percent change & statistics
    fire_rate_mean_morn = mean(fire_rate_morn{xx}, 'omitnan');
    fire_rate_mean_noon = mean(fire_rate_noon{xx}, 'omitnan');
    if strcmp(effect_sz_test, 'Perc_Change')
        fire_rate_effect_size = ((fire_rate_mean_noon - fire_rate_mean_morn) / fire_rate_mean_morn) * 100;
    elseif strcmp(effect_sz_test, 'Cohen')
        cohen_d_effect = meanEffectSize(fire_rate_morn{xx}, fire_rate_noon{xx}, Effect = 'cohen');
        fire_rate_effect_size = cohen_d_effect.Effect;
    end
    fprintf('Effect size is %0.2f \n', fire_rate_effect_size)

    %% Save the file if selected
    Fig_Title = strcat(Fig_Title, {' '}, 'Scatter');
    Save_Figs(Fig_Title, Save_File)

end % End the xds loop







