%% Loading the morning & afternoon files
clear
clc

% Monkey Name
Monkey = 'Pop';
% Select the date & task to analyze (YYYYMMDD)
Date = '20211001';
Task = 'PG';
% Do you want to process the XDS file? (1 = yes; 0 = no)
Process_XDS = 1;

xds_morn = Load_XDS(Monkey, Date, Task, 'Morn', Process_XDS);
xds_noon = Load_XDS(Monkey, Date, Task, 'Noon', Process_XDS);

[xds_morn, xds_noon] = Mismatch(xds_morn, xds_noon, 0);

%% Some variable extraction & definitions

Save_Figs = 0;

unit_names = xds_morn.unit_names;

% Font specifications
label_font_size = 17;
legend_font_size = 13;
title_font_size = 14;
save_title_font_size = 30;
font_name = 'Arial';

% Plotting specifications
marker_metric ='.';
sz = 500;

%% Calculate the mean & std nonlinear energy

avg_nonlin_morn = zeros(length(unit_names),1);
std_nonlin_morn = zeros(length(unit_names),1);
avg_nonlin_noon = zeros(length(unit_names),1);
std_nonlin_noon = zeros(length(unit_names),1);

for ii = 1:length(unit_names)
    avg_nonlin_morn(ii) = mean(xds_morn.nonlin_waveforms{1,ii}, 'All');
    std_nonlin_morn(ii) = std(xds_morn.nonlin_waveforms{1,ii}, 0, 'All');
    avg_nonlin_noon(ii) = mean(xds_noon.nonlin_waveforms{1,ii}, 'All');
    std_nonlin_noon(ii) = std(xds_noon.nonlin_waveforms{1,ii}, 0, 'All');
end

%% Plot the individual nonlinear energy

sort_p_value = zeros(length(unit_names),1);

for ii = 1:length(unit_names)

    [~, sort_p_value(ii)] = NonLinearEnergy_Morn_v_Noon(xds_morn, xds_noon, char(unit_names(ii)), 0, 0);

    % Define the save directory & save the figures
    if ~isequal(Save_Figs, 0)
        save_dir = 'C:\Users\rhpow\Documents\Grad School\Pop\20210922\XDS\Unsorted\PG\Nonlinear Energy Figures\';
        fig_title = strcat('Nonlinear Energy - ', '', unit_names(ii));
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(fig_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'png')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'pdf')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'fig')
        end
        close all
    end

end

%% Plot the total nonlinear energy

figure
hold on

% Set the title

title('Morning vs. Afternoon Nonlinear Energy', 'FontSize', title_font_size)

% Label the axis
xlabel('Morning Nonlinear Energy', 'FontSize', label_font_size);
ylabel('Afternoon Nonlinear Energy', 'FontSize', label_font_size);

for ii = 1:length(unit_names)

    
    if sort_p_value(ii) <= 0.05 % Significant change in nonlinear energy
        color_metric = 0;
    elseif sort_p_value(ii) > 0.05 % Insignificant change in nonlinear energy
        color_metric = 1;
    end

    scatter(avg_nonlin_morn(ii), avg_nonlin_noon(ii), sz, marker_metric, 'MarkerEdgeColor', ... 
                [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);

    % Label the unit
    %text(avg_nonlin_morn(ii) + 1000, avg_nonlin_noon(ii) - 1000, extractAfter(unit_names(ii),"elec"));

end

% Plot dummy points for the right legend
dummy_red = plot(-3200, -3200, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor',[1, 0, 0], 'MarkerFaceColor',[1, 0, 0], 'LineWidth', 1.5);
dummy_black = plot(-4000, -4000, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);

% Calculate the axis limit
min_axis_morn = min(avg_nonlin_morn);
min_axis_noon = min(avg_nonlin_noon);
axis_min = round(min(min_axis_morn, min_axis_noon)/5)*5;
max_axis_morn = max(avg_nonlin_morn);
max_axis_noon = max(avg_nonlin_noon);
axis_max = round(max(max_axis_morn, max_axis_noon)/5)*5;

axis_zoom = 4000;
% Draw the identity line 
line([axis_min - axis_zoom, axis_max + axis_zoom],[axis_min - axis_zoom, axis_max + axis_zoom], ... 
    'Color', 'k', 'Linewidth', 1, 'Linestyle','--')

% Plot the top left legend
left_legend = legend([dummy_red, dummy_black], ... 
    {'p > 0.05', 'p ≤ 0.05'}, ... 
    'FontSize', legend_font_size, 'Location', 'NorthWest');
left_legend.ItemTokenSize(1) = 10;

legend box off

xlim([axis_min - axis_zoom, axis_max + axis_zoom])
ylim([axis_min - axis_zoom, axis_max + axis_zoom])

% Only label every other tick
figure_axes = gca;
figure_axes.XAxis.Exponent = 0;
x_labels = string(figure_axes.XAxis.TickLabels);
x_labels(2:2:end) = NaN;
figure_axes.XAxis.TickLabels = x_labels;
figure_axes.YAxis.Exponent = 0;
y_labels = string(figure_axes.YAxis.TickLabels);
y_labels(2:2:end) = NaN;
figure_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off')
% Set The Font
set(figure_axes,'fontname', font_name);

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


