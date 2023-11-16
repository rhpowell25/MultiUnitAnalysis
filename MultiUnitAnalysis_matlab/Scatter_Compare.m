function Scatter_Compare(Monkey, Sampling_Params, Save_File)

%% Load the output structures

Sampling_Params.trial_task = 'PG';
[xds_depth_excel_first, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel_first, column_names_first] = Split_Depth_Excel(xds_depth_excel_first);

Sampling_Params.trial_task = 'KG';
[xds_depth_excel_second, ~] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel_second, column_names_second] = Split_Depth_Excel(xds_depth_excel_second);

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('Baseline', 'Peak', 'Depth', 'Depth_Change')?
fire_rate_phase = 'Depth_Change';

% Do you want to look at the morning or afternoon ('Morn', 'Noon')
Morn_vs_Noon = 'Morn';

% Do you want the name of each unit labeled? (1 = Yes, 0 = No)
unit_label = 1;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_File, 0)
    close all
end

%% Reassign variables according to what you're plotting

if strcmp(fire_rate_phase, 'Baseline')
    disp('Baseline Firing Rate')
    if strcmp(Morn_vs_Noon, 'Morn')
        fire_rate_first = split_depth_excel_first{strcmp(column_names_first, 'bsfr_morn')};
        fire_rate_second = split_depth_excel_second{strcmp(column_names_second, 'bsfr_morn')};
    elseif strcmp(Morn_vs_Noon, 'Noon')
        fire_rate_first = split_depth_excel_first{strcmp(column_names_first, 'bsfr_noon')};
        fire_rate_second = split_depth_excel_second{strcmp(column_names_second, 'bsfr_noon')};
    end
end
if strcmp(fire_rate_phase, 'Peak')
    disp('Peak Firing Rate')
    if strcmp(Morn_vs_Noon, 'Morn')
        bsfr_first = split_depth_excel_first{strcmp(column_names_first, 'bsfr_morn')};
        bsfr_second = split_depth_excel_second{strcmp(column_names_second, 'bsfr_morn')};
        depth_first = split_depth_excel_first{strcmp(column_names_first, 'depth_morn')};
        depth_second = split_depth_excel_second{strcmp(column_names_second, 'depth_morn')};
        fire_rate_first = struct([]);
        fire_rate_second = struct([]);
        for ii = 1:length(bsfr_first)
            fire_rate_first{ii} = bsfr_first{ii} + depth_first{ii};
            fire_rate_second{ii} = bsfr_second{ii} + depth_second{ii};
        end
    elseif strcmp(Morn_vs_Noon, 'Noon')
        bsfr_first = split_depth_excel_first{strcmp(column_names_first, 'bsfr_noon')};
        bsfr_second = split_depth_excel_second{strcmp(column_names_second, 'bsfr_noon')};
        depth_first = split_depth_excel_first{strcmp(column_names_first, 'depth_noon')};
        depth_second = split_depth_excel_second{strcmp(column_names_second, 'depth_noon')};
        fire_rate_first = struct([]);
        fire_rate_second = struct([]);
        for ii = 1:length(bsfr_first)
            fire_rate_first{ii} = bsfr_first{ii} + depth_first{ii};
            fire_rate_second{ii} = bsfr_second{ii} + depth_second{ii};
        end
    end
end
if strcmp(fire_rate_phase, 'Depth')
    disp('Depth of Modulation')
    if strcmp(Morn_vs_Noon, 'Morn')
        fire_rate_first = split_depth_excel_first{strcmp(column_names_first, 'depth_morn')};
        fire_rate_second = split_depth_excel_second{strcmp(column_names_second, 'depth_morn')};
    elseif strcmp(Morn_vs_Noon, 'Noon')
        fire_rate_first = split_depth_excel_first{strcmp(column_names_first, 'depth_noon')};
        fire_rate_second = split_depth_excel_second{strcmp(column_names_second, 'depth_noon')};
    end
end

if strcmp(fire_rate_phase, 'Depth_Change')
    disp('Depth of Modulation Change')
    depth_first_morn = split_depth_excel_first{strcmp(column_names_first, 'depth_morn')};
    depth_second_morn = split_depth_excel_second{strcmp(column_names_second, 'depth_morn')};
    depth_first_noon = split_depth_excel_first{strcmp(column_names_first, 'depth_noon')};
    depth_second_noon = split_depth_excel_second{strcmp(column_names_second, 'depth_noon')};
    fire_rate_first = struct([]);
    fire_rate_second = struct([]);
    for ii = 1:length(depth_first_morn)
        fire_rate_first{ii} = depth_first_noon{ii} - depth_first_morn{ii};
        fire_rate_second{ii} = depth_second_noon{ii} - depth_second_morn{ii};
    end
end

% Extract the other variables
drug_dose = split_depth_excel_first{strcmp(column_names_first, 'drug_dose_mg_per_kg')};
unit_names_first = split_depth_excel_first{strcmp(column_names_first, 'unit_names')};
unit_names_second = split_depth_excel_second{strcmp(column_names_second, 'unit_names')};

%% Some variable extraction & definitions

% Font specifications
label_font_size = 30;
title_font_size = 14;
fig_size = 700;

% Scatter Marker Shapes
marker_metric = '.';

% Scatter Marker Colors
color_metric = 0;

% Scatter Marker sizes
sz = 500;

%% Loop through each of the experimental sessions
for xx = 1:length(xds_depth_excel_first)

    %% Use only shared units

    all_unit_names = intersect(unit_names_first{xx}, unit_names_second{xx});

    %% Add the monkey name to the title

    if strcmp(Monkey, 'All')
        Fig_Title = strcat('All Monkeys,', {' '});
    else
        Fig_Title = '';
        for ii = 1:length(Monkey)
            Fig_Title = strcat(Fig_Title, Monkey{ii}, ',', {' '});
        end
    end

    %% Plot the Depth of Modulation Scatter

    scatter_fig = figure;
    scatter_fig.Position = [200 50 fig_size fig_size];
    hold on

    % Set the title
    if strcmp(Sampling_Params.trial_sessions, 'Ind')
        title(strcat(Fig_Title, file_names{xx}, {' '}, string(drug_dose{xx}(1)), {', '}, Morn_vs_Noon), 'FontSize', title_font_size)
    elseif strcmp(Sampling_Params.trial_sessions, 'All')
        scatter_title = 'All Trials';
        title(strcat(Fig_Title, scatter_title, {' '}, Drug_Choice), 'FontSize', title_font_size)
    end

    % Label the axis
    if strcmp(fire_rate_phase, 'Baseline')
        xlabel('First Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
        ylabel('Second Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
    end
    if strcmp(fire_rate_phase, 'Peak')
        xlabel('First Peak Firing Rate (Hz)', 'FontSize', label_font_size);
        ylabel('Second Peak Firing Rate (Hz)', 'FontSize', label_font_size);
    end
    if strcmp(fire_rate_phase, 'Depth')
        xlabel('First Depth of Modulation (Hz)', 'FontSize', label_font_size);
        ylabel('Second Depth of Modulation (Hz)', 'FontSize', label_font_size);
    end
    if strcmp(fire_rate_phase, 'Depth_Change')
        xlabel('First Change in Modulation (Hz)', 'FontSize', label_font_size);
        ylabel('Second Change in Modulation (Hz)', 'FontSize', label_font_size);
    end

    % Calculate the axis limits
    min_first = min(fire_rate_first{xx});
    min_second = min(fire_rate_second{xx});
    axis_min = round(min(min_first, min_second)/5)*5;
    max_first = max(fire_rate_first{xx});
    max_second = max(fire_rate_second{xx});
    axis_max = round(max(max_first, max_second)/5)*5;

    % Draw the unity line 
    line([axis_min - 10, axis_max + 10],[axis_min - 10, axis_max + 10], ... 
        'Color', 'k', 'Linewidth', 1, 'Linestyle','--')

    for jj = 1:length(all_unit_names)

        first_unit_idx = find(strcmp(unit_names_first{xx}, cell2mat(all_unit_names(jj))));
        second_unit_idx = find(strcmp(unit_names_second{xx}, cell2mat(all_unit_names(jj))));

        scatter(fire_rate_first{xx}(first_unit_idx), fire_rate_second{xx}(second_unit_idx), sz, marker_metric, 'MarkerEdgeColor', ... 
            [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);
        if isequal(unit_label, 1)
            text(fire_rate_first{xx}(first_unit_idx) + .5, fire_rate_second{xx}(second_unit_idx) - .5, ...
                extractAfter(all_unit_names(jj), "elec"));
        end

    end % End of unit loop

    % Display the means
    fire_rate_mean_first = mean(fire_rate_first{xx});
    fire_rate_mean_second = mean(fire_rate_second{xx});
    fprintf('First average is %0.2f \n', fire_rate_mean_first)
    fprintf('Second average is %0.2f \n', fire_rate_mean_second)

    % Set the axis
    xlim([axis_min - 10, axis_max + 10])
    ylim([axis_min - 10, axis_max + 10])

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end % End the xds loop


