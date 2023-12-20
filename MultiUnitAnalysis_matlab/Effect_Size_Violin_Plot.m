function Effect_Size_Violin_Plot(Monkey, Sampling_Params, Save_File)

%% Load the output structures

[xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params);
[split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel);

%% Some of the plotting specifications

% Which firing rate phase do you want to plot? ('Baseline', 'Ramp', 'TgtHold', 'Depth')?
fire_rate_phase = 'Depth';

% Do you want to remove outliers? (1 = Yes, 0 = No)
remove_outliers = 1;

% Do you want to manually set the y-axis?
man_y_axis = 'No';
%man_y_axis = [-0.7, 0.3];

% What effect size meausure do you want to use ('Perc', 'Cohen')
effect_sz_test = 'Cohen';

% Do you want the simplified title? (1 = Yes, 0 = No)
simp_title = 0;

% Save the figures to your desktop? ('All', 'pdf', 'png', 'fig', 0 = No)
if ~isequal(Save_File, 0)
    close all
end

%% Reassign variables according to what you're plotting

if strcmp(fire_rate_phase, 'Baseline')
    disp('Baseline Firing Rate')
    if strcmp(effect_sz_test, 'Perc')
        effect_sizes = split_depth_excel{strcmp(column_names, 'bsfr_perc')};
    elseif strcmp(effect_sz_test, 'Cohen')
        effect_sizes = split_depth_excel{strcmp(column_names, 'bsfr_cohen_d')};
    end
end
if strcmp(fire_rate_phase, 'Ramp')
    disp('Ramp Phase')
    if strcmp(effect_sz_test, 'Perc')
        effect_sizes = split_depth_excel{strcmp(column_names, 'ramp_perc')};
    elseif strcmp(effect_sz_test, 'Cohen')
        effect_sizes = split_depth_excel{strcmp(column_names, 'ramp_cohen_d')};
    end
end
if strcmp(fire_rate_phase, 'TgtHold')
    disp('TgtHold Phase')
    if strcmp(effect_sz_test, 'Perc')
        effect_sizes = split_depth_excel{strcmp(column_names, 'TgtHold_perc')};
    elseif strcmp(effect_sz_test, 'Cohen')
        effect_sizes = split_depth_excel{strcmp(column_names, 'TgtHold_cohen_d')};
    end
end
if strcmp(fire_rate_phase, 'Depth')
    disp('Depth of Modulation')
    if strcmp(effect_sz_test, 'Perc')
        effect_sizes = split_depth_excel{strcmp(column_names, 'depth_perc')};
    elseif strcmp(effect_sz_test, 'Cohen')
        effect_sizes = split_depth_excel{strcmp(column_names, 'depth_cohen_d')};
    end
end

% Extract the other variables
drug_dose = split_depth_excel{strcmp(column_names, 'drug_dose_mg_per_kg')};
all_unit_names = split_depth_excel{strcmp(column_names, 'unit_names')};

%% Some variable extraction & definitions

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;
effect_size_dims = [0.7 0.375 0.44 0.44];
if isequal(remove_outliers, 1)
    outlier_ann_dims = [0.7 0.325 0.44 0.44];
end

%% Loop through each of the experimental sessions
for xx = 1:length(all_unit_names)

    % Skip the function if no units match
    if isempty(all_unit_names{xx})
        disp('No units match criteria')
        continue
    end

    %% Title information

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

    %% Plot the Violin Plot

    % Remove outliers
    if isequal(remove_outliers, 1)
        effect_outliers = isoutlier(effect_sizes{xx,1}, 'mean');
        %if ~isnan(depth_exclusion)
        %    num_outliers = length(find(effect_outliers == 1)) + length(depth_violations{xx,1});
        %elseif isnan(depth_exclusion)
            num_outliers = length(find(effect_outliers == 1));
        %end
        effect_sizes{xx,1}(effect_outliers) = [];
    end

    % Define the violin plot color
    if strcmp(Sampling_Params.drug_choice, 'Caff') || strcmp(Sampling_Params.drug_choice, 'Lex')
        violin_color = [0 0.5 0];
    elseif strcmp(Sampling_Params.drug_choice, 'Cyp')
        violin_color = [1 0 0];
    elseif strcmp(Sampling_Params.drug_choice, 'Con')
        violin_color = [0 0 0];
    elseif contains(file_names{1,1}, 'Caff') || contains(file_names{1,1}, 'Lex')
        violin_color = [0 0.5 0];
    elseif contains(file_names{1,1}, 'Cyp')
        violin_color = [1 0 0];
    elseif contains(file_names{1,1}, 'Con')
        violin_color = [0 0 0];
    end

    violin_fig = figure;
    violin_fig.Position = [200 50 Plot_Params.fig_size Plot_Params.fig_size];
    hold on
    effect_positions = (1:length(effect_sizes{xx,1}));
    Violin_Plot(effect_sizes(xx,1), effect_positions, 'ViolinColor', violin_color);

    % Set the axis
    if ~ischar(man_y_axis)
        ylim(man_y_axis)
    end
    xlim([0.5, 1.5])

    ylabel('Effect Size', 'FontSize', Plot_Params.label_font_size)

    % Draw the unity line 
    line([0.5, 1.5], [0, 0], 'Color', 'k', 'Linewidth', Plot_Params.mean_line_width, 'Linestyle','--')

    % Set the title
    title(Fig_Title, 'FontSize', Plot_Params.title_font_size, 'Color', Plot_Params.title_color, 'Interpreter', 'none')
    title('')

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

    % Only label every other tick
    x_labels = string(figure_axes.XAxis.TickLabels);
    y_labels = string(figure_axes.YAxis.TickLabels);
    x_labels(1:end) = NaN;
    y_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    figure_axes.YAxis.TickLabels = y_labels;

    % Find the mean firing rates
    avg_effect_size = mean(effect_sizes{xx,1},'omitnan');
    if strcmp(effect_sz_test, 'Perc')
        effect_size_string = strcat('Δ% =', {' '}, mat2str(round(avg_effect_size*100, 1)), '%');
    elseif strcmp(effect_sz_test, 'Cohen')
        effect_size_string = strcat('Cd =', {' '}, mat2str(round(avg_effect_size, 2)));
    end

    % Annotation of the effect size
    effect_ann_string = {char(effect_size_string)};
    effect_ann = annotation('textbox', effect_size_dims, 'String', effect_ann_string, ... 
        'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
        'EdgeColor','none', 'horizontalalignment', 'Left');
    effect_ann.FontSize = Plot_Params.legend_size;
    effect_ann.FontName = Plot_Params.font_name;

    % Annotation of the outliers removed
    if isequal(remove_outliers, 1)
        if ~ischar(man_y_axis)
            num_outliers = num_outliers + length(find(effect_sizes{xx,1} > man_y_axis(2))) + ...
                length(find(effect_sizes{xx,1} < man_y_axis(1)));
        end
        outlier_string = strcat(mat2str(num_outliers), {' '}, '> 3σ');
        outlier_ann_string = {char(outlier_string)};
        outlier_ann = annotation('textbox', outlier_ann_dims, 'String', outlier_ann_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'Left');
        outlier_ann.FontSize = Plot_Params.legend_size;
        outlier_ann.FontName = Plot_Params.font_name;
    end

    %% Save the file if selected
    Fig_Title = strcat(Fig_Title, {' '}, 'Violin');
    Save_Figs(Fig_Title, Save_File)
    
end % End the xds loop







