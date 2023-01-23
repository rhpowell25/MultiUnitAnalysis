%% Load the output structures
clearvars -except xds_morn & xds_noon & xds & unit_name
clc

% What Drug Do You Want To Load? ('Caff', 'Cyp', 'Lex', 'Con')
Drug_Choice = 'Con';

[~, ~, depth_morn, depth_noon, ~, ~, ~, ~, depth_p_value, ~, nonlin_p_value, fract_contam, ...
    ~, ~, file_name, all_unit_names] = Fast_Load_Depth(Drug_Choice);

% Save the matrix to Excel? (1 = Yes, 0 = No)
Save_Excel = 0;

% Create the matrix headers
sort_stats = cell(length(all_unit_names),7);
sort_stats{1,1} = 'File';
sort_stats{1,2} = 'All Units';
sort_stats{1,3} = '# Fract Contam < 10%';
sort_stats{1,4} = '# Insignificant Nonlinear Energy';
sort_stats{1,5} = '# Both';
sort_stats{1,6} = '# Significant Positive Depth Change';
sort_stats{1,7} = '# Significant Negative Depth Change';

%% Loop through each of the experimental sessions
for xx = 1:length(all_unit_names)

    % Changes in depth of modulation
    depth_change = depth_noon{xx,1} - depth_morn{xx,1};
    pos_depth_change_idx = find(depth_change > 0);
    neg_depth_change_idx = find(depth_change < 0);

    sort_stats{xx + 1,1} = file_name{xx,1};
    sort_stats{xx + 1,2} = length(all_unit_names{xx,1});

    % Find the indexes of the well isolated units
    good_fract_contam_idx = find(fract_contam{xx,1} < 0.1);
    sort_stats{xx + 1,3} = length(good_fract_contam_idx);

    % Find the indexes of the consistently sorted units
    insig_nonlin_p_value_idx = find(nonlin_p_value{xx,1} >= 0.05);
    sort_stats{xx + 1,4} = length(insig_nonlin_p_value_idx);

    % Find the indexes of the consistently well isolated units
    well_sorted_idx = intersect(good_fract_contam_idx, insig_nonlin_p_value_idx);
    sort_stats{xx + 1,5} = length(well_sorted_idx);

    % Find the units whose depth of modulation significantly changed
    sig_depth_p_value_idx = find(depth_p_value{xx,1} < 0.05);
    pos_sig_depth_idx = intersect(pos_depth_change_idx, sig_depth_p_value_idx);
    neg_sig_depth_idx = intersect(neg_depth_change_idx, sig_depth_p_value_idx);
    sort_stats{xx + 1,6} = length(pos_sig_depth_idx);
    sort_stats{xx + 1,7} = length(neg_sig_depth_idx);
    
end

%% Save to Excel

if isequal(Save_Excel, 1)

    % Save the file
    writecell(sort_stats, strcat(Save_Path, strcat('Sort_Stats_', Drug_Choice, '.xlsx')))

end




