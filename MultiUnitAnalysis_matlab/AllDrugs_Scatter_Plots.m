%% Load the output structures
clear
clc

% What Monkey Do You Want To Load? ('All', 'Mihili', 'Jango', 'Jaco', 'Pop', 'Pancake')
Monkey_Choice = strings;
Monkey_Choice{1} = 'Pop';
%Monkey_Choice{1} = 'Tot';
%Monkey_Choice{3} = 'Groot';
%Monkey_Choice{4} = 'Pancake';
% Sorted or unsorted (1 vs 0)
Sorted = 1;

event = 'window_trial_gocue';

Sampling_Params = struct( ...
    'trial_task', 'PG', ... % Which task do you want ('PG', 'WS', 'KG', 'All')?
    'trial_sessions', 'All', ... % Individual Sessions or All Sessions? ('Ind' vs 'All')
    'unit_quality', 'Stable', ... % Unit quality ('All' vs. 'Stable')
    'ISI_quality', 'All', ... % ISI quality ('All' vs. 'Single')
    'depth_change', NaN, ... % What change in depth of modulation do you want to observe (# vs NaN)
    'pref_dir', 'All', ... % What preferred direction do you want to plot(-90, 0, 90, 180, 'All')
    'depth_min', NaN, ... % What minimum depth of modulation do you want to observe (# vs NaN)
    'depth_max', NaN, ... % What maximum depth of modulation do you want to observe (# vs NaN)
    'peaktopeak_min', NaN, ... % What minimum peak to peak amplitude do you want to observe (# vs NaN)
    'mod_sig', 'All', ... % Unit modulation significance ('All' vs. 'Sig')
    'depth_sig', 'All', ... % Unit depth change significance ('All' vs. 'Sig')
    'spike_width', 'All', ... % What spike width ('All', 'Small', 'Medium', 'Large')
    'post_spike_facil', 'All', ... % What level of post-spike facilitation ('All', 'Weak', 'Moderate', 'Strong')
    'max_vs_min', 0); % Remove sessions with only a single target? (1 = Yes, 0 = No)

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('Baseline', 'Ramp', 'TgtHold', 'Depth')?
fire_rate_phase = 'Depth';

% Plot each unit or the averages only? (1 = each unit; 0 = averages only)
each_unit = 0;

% Drug Legend? (1 = Yes, 0 = No)
drug_leg = 0;

% Save the figures to your desktop? (1 = Yes, 0 = No)
Save_Figs = 0;

%% Some variable extraction & definitions

% Font Sizes
label_font_size = 17;
legend_font_size = 13;
title_font_size = 24;
save_title_font_size = 24;
font_name = 'Arial';

%% Generate the figure

figure
hold on

% Set the title
title('Depth of Modulation:', 'FontSize', title_font_size)

% Label the axes
xlabel('Morning Depth of Modulation (Hz)', 'FontSize', label_font_size);
ylabel('Afternoon Depth of Modulation (Hz)', 'FontSize', label_font_size);

% Set up a struct for the elliptical error information
ellipsis_specs = zeros(4,4);

for dd = 1:4
    if dd == 1
        %drug_color = [0.5, 0.5, 1]; % Blue
        drug_color = 'r';
        Drug_Choice = 'Cyp';
        [xds_depth_excel, file_names] = Load_Depth_Excel(Drug_Choice, Monkey_Choice, event, Sorted, Sampling_Params);
        [split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);
    end
    if dd == 2
        drug_color = [0 0.5 0]; % Dark Green
        Drug_Choice = 'Caff';
        [xds_depth_excel, file_names] = Load_Depth_Excel(Drug_Choice, Monkey_Choice, event, Sorted, Sampling_Params);
        [split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);
    end
    if dd == 3
        drug_color = [0 0.5 0]; % Dark Green
        Drug_Choice = 'Lex';
        [xds_depth_excel, file_names] = Load_Depth_Excel(Drug_Choice, Monkey_Choice, event, Sorted, Sampling_Params);
        [split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);
    end
    if dd == 4
        drug_color = [0, 0, 0]; % Black
        Drug_Choice = 'Con';
        [xds_depth_excel, file_names] = Load_Depth_Excel(Drug_Choice, Monkey_Choice, event, Sorted, Sampling_Params);
        [split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);
    end

    %% Reassign variables according to what you're plotting

    if strcmp(fire_rate_phase, 'Baseline')
        disp('Baseline Firing Rate')
        fire_rate_morn = split_depth_excel{strcmp(column_names, 'bsfr_morn')};
        fire_rate_noon = split_depth_excel{strcmp(column_names, 'bsfr_noon')};
        fire_rate_err_morn = split_depth_excel{strcmp(column_names, 'bsfr_err_morn')};
        fire_rate_err_noon = split_depth_excel{strcmp(column_names, 'bsfr_err_noon')};
    end
    if strcmp(fire_rate_phase, 'Ramp')
        disp('Ramp Phase')
        fire_rate_morn = split_depth_excel{strcmp(column_names, 'ramp_morn')};
        fire_rate_noon = split_depth_excel{strcmp(column_names, 'ramp_noon')};
        fire_rate_err_morn = split_depth_excel{strcmp(column_names, 'ramp_err_morn')};
        fire_rate_err_noon = split_depth_excel{strcmp(column_names, 'ramp_err_noon')};
    end
    if strcmp(fire_rate_phase, 'TgtHold')
        disp('TgtHold Phase')
        fire_rate_morn = split_depth_excel{strcmp(column_names, 'TgtHold_morn')};
        fire_rate_noon = split_depth_excel{strcmp(column_names, 'TgtHold_noon')};
        fire_rate_err_morn = split_depth_excel{strcmp(column_names, 'TgtHold_err_morn')};
        fire_rate_err_noon = split_depth_excel{strcmp(column_names, 'TgtHold_err_noon')};
    end
    if strcmp(fire_rate_phase, 'Depth')
        disp('Depth of Modulation')
        fire_rate_morn = split_depth_excel{strcmp(column_names, 'depth_morn')};
        fire_rate_noon = split_depth_excel{strcmp(column_names, 'depth_noon')};
        fire_rate_err_morn = split_depth_excel{strcmp(column_names, 'mp_err_morn')};
        fire_rate_err_noon = split_depth_excel{strcmp(column_names, 'mp_err_noon')};
    end

    %% Plot the elliptical error probable

    % The percent of points you want
    err_percent = .5;
    X = fire_rate_morn{1,1}';
    Y = fire_rate_noon{1,1}';

    [meanX, meanY, ~, ~, ~, X_ellipse, Y_ellipse] = ... 
        Ellip_Err_Prob(X, Y, err_percent);

    for ii = 1:length(err_percent)
        plot(X_ellipse(ii,:),Y_ellipse(ii,:), 'Color', drug_color,'LineWidth',2);
    end
    % Put the ellipses info into the output variable
    ellipsis_specs(1,dd) = min(X_ellipse(1,:));
    ellipsis_specs(2,dd) = max(X_ellipse(1,:));
    ellipsis_specs(3,dd) = min(Y_ellipse(1,:));
    ellipsis_specs(4,dd) = max(Y_ellipse(1,:));

    % Plot the mean
    plot(meanX,meanY, '*', 'MarkerSize',10, ...
            'MarkerEdgeColor', drug_color, 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);

    %% Plot the Depth of Modulation Scatter

    if isequal(each_unit, 1)
        for jj = 1:total_units

            % If the ISI's mode is < 5
            if round(all_trials_merged_ISI_mode(jj,1)) <= 5
                marker_metric = '*';
                sz = 200;
            else
                marker_metric ='.';
                sz = 500;
            end
    
            % Plot each unit
            scatter(all_trials_depth_morn(jj), all_trials_depth_noon(jj), sz, marker_metric, 'MarkerEdgeColor', ... 
                drug_color, 'MarkerFaceColor', drug_color, 'LineWidth', 1.5);

            % Error
            err_morn = errorbar(all_trials_depth_morn(jj), all_trials_depth_noon(jj), ... 
                all_trials_mpfr_err_morn(jj), 'horizontal');
            err_morn.Color = drug_color;
            err_morn.LineWidth = 1;
            err_noon = errorbar(all_trials_depth_morn(jj), all_trials_depth_noon(jj), ... 
                all_trials_mpfr_err_noon(jj), 'vertical');
            err_noon.Color = drug_color;
            err_noon.LineWidth = 1;
        end

    end % End of unit loop
    
end % End of drug loop

% Set the axis
if isequal(each_unit, 1)
    min_err_morn = min(all_trials_depth_morn - all_trials_mpfr_err_morn);
    min_err_noon = min(all_trials_depth_noon - all_trials_mpfr_err_noon);
    axis_min = round(min(min_err_morn, min_err_noon)/5)*5;
    max_err_morn = max(all_trials_depth_morn + all_trials_mpfr_err_morn);
    max_err_noon = max(all_trials_depth_noon + all_trials_mpfr_err_noon);
    axis_max = round(max(max_err_morn, max_err_noon)/5)*5;
else
    axis_min = min(ellipsis_specs, [], 'All');
    axis_max = max(ellipsis_specs, [], 'All');
end

xlim([axis_min - 5, axis_max + 5])
ylim([axis_min - 5, axis_max + 5])

% Draw the identity line 
line([axis_min - 5, axis_max + 5],[axis_min - 5, axis_max + 5], ... 
    'Color', 'k', 'Linewidth', 1, 'Linestyle','--')

if isequal(drug_leg, 1)
    % Plot dummy points for the legend
    dummy_orange = plot(-12,-12, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'MarkerFaceColor',[1, 0, 0], 'LineWidth', 1.5);
    dummy_red = plot(-17,-17, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', [0.6350 0.0780 0.1840], 'MarkerFaceColor',[0.5, 0, 0], 'LineWidth', 1.5);
    dummy_blue = plot(-20,-20, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor',[0.5, 0.5, 1], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);
    dummy_black = plot(-22,-22, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);
        
    % Plot the legend
    legend([dummy_orange, dummy_red, dummy_blue, dummy_black], ... 
        {'Lexapro', 'Caffeine', 'Cyproheptadine', 'Control'}, ... 
        'FontSize', legend_font_size, 'Location', 'southeast')
    legend boxoff

end

% Only label every other tick
figure_axes = gca;
x_labels = string(figure_axes.XAxis.TickLabels);
y_labels = string(figure_axes.YAxis.TickLabels);
x_labels(2:2:end) = NaN;
y_labels(2:2:end) = NaN;
figure_axes.XAxis.TickLabels = x_labels;
figure_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off')
% Set The Font
set(figure_axes,'fontname', font_name);

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
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
        title(fig_title, 'Fontsize', save_title_font_size);
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, 'vs.', 'vs');
        save_title = strrep(save_title, 'mg.', 'mg');
        save_title = strrep(save_title, 'kg.', 'kg');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
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

            




