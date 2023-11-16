function Scatter_PrefDir(Monkey, Sampling_Params, Save_File)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('Baseline', 'Peak', 'Depth', 'Depth_Change')?
fire_rate_phase = 'Depth';

% Do you want to look at the morning or afternoon ('Morn', 'Noon')
Morn_vs_Noon = 'Morn';

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_File, 0)
    close all
end

%% Reassign variables according to what you're plotting
if strcmp(Morn_vs_Noon, 'Morn')
    avg_bsfr = split_depth_excel{strcmp(column_names, 'bsfr_morn')};
elseif strcmp(Morn_vs_Noon, 'Noon')
    avg_bsfr = split_depth_excel{strcmp(column_names, 'bsfr_noon')};
end

if strcmp(fire_rate_phase, 'Baseline')
    disp('Baseline Firing Rate')
    if strcmp(Morn_vs_Noon, 'Morn')
        fire_rate_prefdir = split_depth_excel{strcmp(column_names, 'bsfr_morn')};
    elseif strcmp(Morn_vs_Noon, 'Noon')
        fire_rate_prefdir = split_depth_excel{strcmp(column_names, 'bsfr_noon')};
    end
end
if strcmp(fire_rate_phase, 'Peak')
    disp('Peak Firing Rate')
    if strcmp(Morn_vs_Noon, 'Morn')
        depth_prefdir = split_depth_excel{strcmp(column_names, 'depth_morn')};
        fire_rate_prefdir = struct([]);
        for ii = 1:length(avg_bsfr)
            fire_rate_prefdir{ii} = avg_bsfr{ii} + depth_prefdir{ii};
        end
    elseif strcmp(Morn_vs_Noon, 'Noon')
        depth_prefdir = split_depth_excel{strcmp(column_names, 'depth_noon')};
        fire_rate_prefdir = struct([]);
        for ii = 1:length(avg_bsfr)
            fire_rate_prefdir{ii} = avg_bsfr{ii} + depth_prefdir{ii};
        end
    end
end
if strcmp(fire_rate_phase, 'Depth')
    disp('Depth of Modulation')
    if strcmp(Morn_vs_Noon, 'Morn')
        fire_rate_prefdir = split_depth_excel{strcmp(column_names, 'depth_morn')};
    elseif strcmp(Morn_vs_Noon, 'Noon')
        fire_rate_prefdir = split_depth_excel{strcmp(column_names, 'depth_noon')};
    end
end

if strcmp(fire_rate_phase, 'Depth_Change')
    disp('Depth of Modulation Change')
    depth_prefdir_morn = split_depth_excel{strcmp(column_names, 'depth_morn')};
    depth_prefdir_noon = split_depth_excel{strcmp(column_names, 'depth_noon')};
    fire_rate_prefdir = struct([]);
    for ii = 1:length(depth_prefdir_morn)
        fire_rate_prefdir{ii} = depth_prefdir_noon{ii} - depth_prefdir_morn{ii};
    end
end

% Extract the other variables
unit_names = split_depth_excel{strcmp(column_names, 'unit_names')};
pref_dirs = split_depth_excel{strcmp(column_names, 'pref_dir')};
targets = split_depth_excel{strcmp(column_names, 'target')};
if strcmp(Morn_vs_Noon, 'Morn')
    max_fr_time = split_depth_excel{strcmp(column_names, 'alignment_morn')};
elseif strcmp(Morn_vs_Noon, 'Noon')
    max_fr_time = split_depth_excel{strcmp(column_names, 'alignment_noon')};
end

%% Extract the firing rate from the opposite preferred direction

fire_rate_oppdir = struct([]);
for ii = 1:length(file_names)
    % load the corresponding xds file

    xtra_info = extractAfter(file_names{ii}, ' ');
    % Date
    Date = erase(file_names{ii}, strcat({' '}, xtra_info));

    % Load the xds files
    xds = Load_XDS(Monkey, Date, Sampling_Params.trial_task, ...
        Sampling_Params.sorted, Morn_vs_Noon);

    % Pull the binning paramaters
    [Bin_Params] = Binning_Parameters;
    % Window to calculate max firing rate
    half_window_length = Bin_Params.half_window_length; % Time (sec.)

    % Loop through each unit
    for pp = 1:length(unit_names{ii})
        unit_name = unit_names{ii}{pp};
        pref_dir = pref_dirs{ii}(pp);
        target = targets{ii}(pp);

        % Define the opposite of the preferred direction
        if isequal(abs(pref_dir), 90)
            opp_dir = pref_dir*-1;
        elseif isequal(pref_dir, 0)
            opp_dir = 180;
        elseif isequal(pref_dir, 180)
            opp_dir = 0;
        end

        % Times for rewarded trials
        [Alignment_Times] = EventAlignmentTimes(xds, opp_dir, target, Sampling_Params.event);
        % Find the unit of interest
        [N] = Find_Unit(xds, unit_name);
        % Extract all the spikes of the unit
        spikes = xds.spikes{1, N};

        % Get the movement phase firing rates
        mp_fr = zeros(length(Alignment_Times),1);
        for tt = 1:length(mp_fr)
            t_start = Alignment_Times(tt) + max_fr_time{ii}(pp) - half_window_length;
            t_end = Alignment_Times(tt) + max_fr_time{ii}(pp) + half_window_length;
            mp_fr(tt,1) = length(find((spikes >= t_start) & ...
                    (spikes <= t_end))) / (2*half_window_length);
        end

        % Find the depth of modulation
        fire_rate_oppdir{ii}(pp,1) = mean(mp_fr) - avg_bsfr{ii}(pp,1);

    end

end

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
for xx = 1:length(xds_depth_excel)

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
        title(strcat(Fig_Title, file_names{xx}, {' '}, Morn_vs_Noon), 'FontSize', title_font_size)
    elseif strcmp(Sampling_Params.trial_sessions, 'All')
        scatter_title = 'All Trials';
        title(strcat(Fig_Title, scatter_title, {' '}, Drug_Choice), 'FontSize', title_font_size)
    end

    % Label the axis
    if strcmp(fire_rate_phase, 'Baseline')
        xlabel('Preferred Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
        ylabel('Opposite Baseline Firing Rate (Hz)', 'FontSize', label_font_size);
    end
    if strcmp(fire_rate_phase, 'Peak')
        xlabel('Preferred Peak Firing Rate (Hz)', 'FontSize', label_font_size);
        ylabel('Opposite Peak Firing Rate (Hz)', 'FontSize', label_font_size);
    end
    if strcmp(fire_rate_phase, 'Depth')
        xlabel('Preferred Depth of Modulation (Hz)', 'FontSize', label_font_size);
        ylabel('Opposite Depth of Modulation (Hz)', 'FontSize', label_font_size);
    end
    if strcmp(fire_rate_phase, 'Depth_Change')
        xlabel('Preferred Change in Modulation (Hz)', 'FontSize', label_font_size);
        ylabel('Opposite Change in Modulation (Hz)', 'FontSize', label_font_size);
    end

    % Calculate the axis limits
    min_prefdir = min(fire_rate_prefdir{xx});
    min_oppdir = min(fire_rate_oppdir{xx});
    axis_min = round(min(min_prefdir, min_oppdir)/5)*5;
    max_prefdir = max(fire_rate_prefdir{xx});
    max_oppdir = max(fire_rate_oppdir{xx});
    axis_max = round(max(max_prefdir, max_oppdir)/5)*5;

    % Draw the unity & zero line 
    line([axis_min - 10, axis_max + 10],[axis_min - 10, axis_max + 10], ... 
        'Color', 'k', 'Linewidth', 1, 'Linestyle','--')
    line([axis_min - 10, axis_max + 10],[0, 0], ... 
        'Color', 'k', 'Linewidth', 1, 'Linestyle','--')

    for jj = 1:length(unit_names)

        scatter(fire_rate_prefdir{xx}, fire_rate_oppdir{xx}, sz, marker_metric, 'MarkerEdgeColor', ... 
            [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);

    end % End of unit loop

    % Display the means
    fire_rate_mean_prefdir = mean(fire_rate_prefdir{xx});
    fire_rate_mean_oppdir = mean(fire_rate_oppdir{xx});
    fprintf('Preferred average is %0.2f \n', fire_rate_mean_prefdir)
    fprintf('Opposite average is %0.2f \n', fire_rate_mean_oppdir)

    % Set the axis
    xlim([axis_min - 10, axis_max + 10])
    ylim([axis_min - 10, axis_max + 10])

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end % End the xds loop


