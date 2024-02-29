 function [morn_and_noon] = Multi_Session_Matrix(Matrix_Length)

%% Create matrix for all units morning and afternoon

% Define the excel matrix headers
matrix_headers = {'unit_names', 'bsfr_morn', 'bsfr_noon', ...
    'depth_morn', 'depth_noon', 'bsfr_err_morn', 'bsfr_err_noon', 'depth_err_morn', ...
    'depth_err_noon', 'bsfr_p_val', 'bsfr_perc', 'bsfr_cohen_d', 'mod_p_val_morn', ...
    'mod_p_val_noon', 'depth_p_val', 'depth_perc', 'depth_cohen_d', 'pref_dir', ...
    'target', 'alignment', 'post_spike_facil', 'wave_sigtonoise', 'wave_peaktopeak', ...
    'wave_p_val', 'nonlin_sigtonoise', 'nonlin_peaktopeak', 'nonlin_p_val', ...
    'spike_width', 'repol_time', 'fract_contam', 'num_targets', 'rxn_time', 'drug_dose_mg_per_kg'};

% Create the excel matrix
morn_and_noon = array2table(zeros(Matrix_Length, length(matrix_headers)));

morn_and_noon.Properties.VariableNames = matrix_headers;




