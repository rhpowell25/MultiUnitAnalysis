function Multi_Session_Depth(Monkey, event, Drug_Choice, Save_Excel)

%% Display the function being used
disp('Multi-Experiment Depth Results:');

%% Some of the analysis specifications

% What minimum signal to noise ration? (# or NaN)
SNR_Minimum = NaN; %5;

% Which targets do you want the mnovement phase firing rate calculated from? ('Max', 'Min', 'All')
tgt_mpfr = 'Max';

% Sorted or unsorted (1 vs 0)
Sorted = 1;

if isequal(Sorted, 1)
    Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '\Sorted\');
else
    Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '\Unsorted\');
end

% Load the file information
[Dates, Tasks, Drug_Dose] = File_Details(Monkey, Drug_Choice);

%% Build the output arrays

% Morning baseline firing rate
bsfr_morn = struct([]);
% Afternoon baseline firing rate
bsfr_noon = struct([]);

% Morning depth
depth_morn = struct([]);
% Afternoon depth
depth_noon = struct([]);
    
% Morning error
bsfr_err_morn = struct([]);
depth_err_morn = struct([]);
% Afternoon error
bsfr_err_noon = struct([]);
depth_err_noon = struct([]);

% Modulation statistics
mod_p_value_t_test_morn = struct([]);
mod_p_value_wilcoxon_morn = struct([]);
mod_p_value_t_test_noon = struct([]);
mod_p_value_wilcoxon_noon = struct([]);
% Baseline firing rate statistics
bsfr_p_value_ks_test_morn = struct([]);
bsfr_p_value_ks_test_noon = struct([]);
bsfr_p_value_t_test = struct([]);
bsfr_p_value_wilcoxon = struct([]);
bsfr_effect_perc = struct([]);
bsfr_effect_cohen_d = struct([]);
% Depth of modulation statistics
depth_p_value_ks_test_morn = struct([]);
depth_p_value_ks_test_noon = struct([]);
depth_p_value_t_test = struct([]);
depth_p_value_wilcoxon = struct([]);
depth_effect_perc = struct([]);
depth_effect_cohen_d = struct([]);

% Signal to noise ratio
wave_sigtonoise = struct([]);
nonlin_sigtonoise = struct([]);
% Peak to peak amplitude
wave_peaktopeak = struct([]);
nonlin_peaktopeak = struct([]);
% Wave shape t-tests
wave_p_value = struct([]);
nonlin_p_value = struct([]);
% Spike Width
spike_width = struct([]);
% Repolarization Time
repol_time = struct([]);
% Fractional Contamination
fract_contam = struct([]);

% Preferred Direction
all_units_pref_dir = struct([]);
% Target center
all_units_target = struct([]);
% Alignment times
alignment = struct([]);

% Post Spike Facilitation
post_spike_facil = struct([]);

% Reaction Time
rxn_time = struct([]);

% Unit Names
all_unit_names = struct([]);

% Number of targets
num_targets = struct([]);

%% Loop through the different experiments
for xx = 1:length(Dates)
    
    % Load the relevant xds file
    xds_morn = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Noon');

    % Process the xds files
    Match_The_Targets = 1;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);

    %[xds_morn] = Subsample_Reaction_Outliers(xds_morn);
    %[xds_noon] = Subsample_Reaction_Outliers(xds_noon);

    % Set the succesful unit counter
    cc = 0;

    %% Find the mean reaction time

    % Which task metric? 'Rxn_Time', 'Trial_Length'
    Task_Metric = 'Rxn_Time';

    [rxn_times] = Task_Metric_ViolinPlot(xds_morn, xds_noon, Task_Metric, 0, 0);

    if length(rxn_times) > 1
        rxn_time{xx,1} = mean(rxn_times);
    else
        rxn_time{xx,1} = rxn_times;
    end

    %% Loop through the units  
    for jj = 1:length(xds_morn.unit_names)

        %% Make sure the morning and noon units are identical
        if strcmp(char(xds_morn.unit_names(jj)), char(xds_noon.unit_names(jj)))
            unit_name = char(xds_morn.unit_names(jj));
            fprintf('%s \n', unit_name);
        else
            disp('Units Are Mismatched')
            return
        end

        %% Check if the unit is stable
        
        [wave_peak_to_peak, wave_spike_width, wave_repol_time, wave_sort_metric] = ...
            WaveShapes_Morn_v_Noon(xds_morn, xds_noon, unit_name, 0, 0);
        if isnan(wave_sort_metric)
            % If the unit doesn't exist or has less than 1000 spikes
            continue
        end
        [nonlin_peak_to_peak, ~, nonlin_sort_metric] = ...
            NonLinearEnergy_Morn_v_Noon(xds_morn, xds_noon, unit_name, 0, 0);

        %% Check if the unit's SNR is below the defined minimum
        
        [wave_SNR_matrix_morn] = SignalToNoise(xds_morn, unit_name, 'Wave');
        wave_SNR_morn = wave_SNR_matrix_morn{2,1};
        [wave_SNR_matrix_noon] = SignalToNoise(xds_noon, unit_name, 'Wave');
        wave_SNR_noon = wave_SNR_matrix_noon{2,1};
        wave_sig_to_noise = (wave_SNR_morn + wave_SNR_noon)/2;
        if ~isnan(SNR_Minimum)
            if wave_SNR_morn <= SNR_Minimum || wave_SNR_noon <= SNR_Minimum
                continue
            end
        end

        [nonlin_SNR_matrix_morn] = SignalToNoise(xds_morn, unit_name, 'Nonlin');
        nonlin_SNR_morn = nonlin_SNR_matrix_morn{2,1};
        [nonlin_SNR_matrix_noon] = SignalToNoise(xds_noon, unit_name, 'Nonlin');
        nonlin_SNR_noon = nonlin_SNR_matrix_noon{2,1};
        nonlin_sig_to_noise = (nonlin_SNR_morn + nonlin_SNR_noon)/2;

        %% Check if the unit is single or multi

        [Fractional_Contam, ~] = InterstimulusInterval_Morn_v_Noon(xds_morn, xds_noon, unit_name, 0, 0);

        [fract_contam_morn, ~] = InterstimulusInterval(xds_morn, unit_name, 0, 0);
        [fract_contam_noon, ~] = InterstimulusInterval(xds_noon, unit_name, 0, 0);

        % Skip the unit if the fract_contam changes between sessions
        if abs(fract_contam_noon - fract_contam_morn) >= 0.5
            continue
        end

        %% Get the baseline firing rates
        [avg_bsfr_morn, ~, err_bsfr_morn, pertrial_bsfr_morn] = ... 
            BaselineFiringRate(xds_morn, unit_name);
        [avg_bsfr_noon, ~, err_bsfr_noon, pertrial_bsfr_noon] = ... 
            BaselineFiringRate(xds_noon, unit_name);

        %% Extract the target directions & centers

        [target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
        [target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);
        
        %% Check to see if both sessions use a consistent number of targets

        % Find matching targets between the two sessions
        [Matching_Idxs_Morn, Matching_Idxs_Noon] = ...
            Match_Targets(target_dirs_morn, target_dirs_noon, target_centers_morn, target_centers_noon);
        
        % Only use the info of target centers conserved between morn & noon
        if ~all(Matching_Idxs_Morn) || ~all(Matching_Idxs_Noon)
            disp('Uneven Targets Between Morning & Afternoon');
            target_centers_morn = target_centers_morn(Matching_Idxs_Morn);
            target_centers_noon = target_centers_noon(Matching_Idxs_Noon);
            target_dirs_morn = target_dirs_morn(Matching_Idxs_Morn);
            target_dirs_noon = target_dirs_noon(Matching_Idxs_Noon);
        end

        %% Only look at the preferred direction
        [pref_dir] = Pref_Direction_Morn_v_Noon(xds_morn, xds_noon, unit_name, event, tgt_mpfr);

        if strcmp(tgt_mpfr, 'Max')
            pref_dir_tgt_morn = max(target_centers_morn(target_dirs_morn == pref_dir));
            pref_dir_tgt_noon = max(target_centers_noon(target_dirs_noon == pref_dir));
        elseif strcmp(tgt_mpfr, 'Min')
            pref_dir_tgt_morn = min(target_centers_morn(target_dirs_morn == pref_dir));
            pref_dir_tgt_noon = min(target_centers_noon(target_dirs_noon == pref_dir));
        end

        if isempty(pref_dir_tgt_morn) || isempty(pref_dir_tgt_noon)
            disp('No targets in the preferred direction!')
            continue
        end

        if isequal(pref_dir_tgt_noon, pref_dir_tgt_noon)
            target_center = pref_dir_tgt_noon;
        end

        %% Find the post-spike facilitation

        %[~, peak_to_noise_ratio_morn] = Spike_Trigger_Avg(xds_morn, unit_name, pref_dir, 0, 0);
        %[~, peak_to_noise_ratio_noon] = Spike_Trigger_Avg(xds_noon, unit_name, pref_dir, 0, 0);

        %max_facil = max(cat(1, peak_to_noise_ratio_morn, peak_to_noise_ratio_noon));
        max_facil = NaN;

        %% Find the number of targets in the preferred direction
        num_targets_morn = length(target_centers_morn(target_dirs_morn == pref_dir));
        num_targets_noon = length(target_centers_morn(target_dirs_morn == pref_dir));
        % Confirm the # of targets in the morning & afternoon are equal
        if isequal(num_targets_morn, num_targets_noon)
            pref_dir_targets = num_targets_morn;
        else
            pref_dir_targets = NaN;
            disp('Unequal number of targets!')
        end

        %% Find the alignment times

        [pertrial_mpfr_morn, pertrial_mpfr_noon, max_fr_time] = ...
            EventWindow_Morn_v_Noon(xds_morn, xds_noon, unit_name, pref_dir, target_center, event);
        alignment_time = max_fr_time;
        % Calculate the depth of modulation per trial
        pertrial_depth_morn = struct([]);
        pertrial_depth_noon = struct([]);
        pertrial_depth_morn{1,1} = pertrial_mpfr_morn{1,1} - avg_bsfr_morn(1,1);
        pertrial_depth_noon{1,1} = pertrial_mpfr_noon{1,1} - avg_bsfr_noon(1,1);
        % Find the mean & standard error of the depth of modulation
        avg_depth_morn = mean(pertrial_depth_morn{1,1});
        avg_depth_noon = mean(pertrial_depth_noon{1,1});
        err_depth_morn = std(pertrial_depth_morn{1,1}) / ...
            sqrt(length(pertrial_depth_morn{1,1}));
        err_depth_noon = std(pertrial_depth_noon{1,1}) / ...
            sqrt(length(pertrial_depth_noon{1,1}));

        %% Check if the unit's depth of modulation changed significantly

        % Modulation statistics & effect sizes
        [~, mod_t_test_morn] = ttest2(pertrial_bsfr_morn{1,1}, pertrial_mpfr_morn{1,1});
        [mod_wilcoxon_morn, ~] = ranksum(pertrial_bsfr_morn{1,1}, pertrial_mpfr_morn{1,1});
        [~, mod_t_test_noon] = ttest2(pertrial_bsfr_noon{1,1}, pertrial_mpfr_noon{1,1});
        [mod_wilcoxon_noon, ~] = ranksum(pertrial_bsfr_noon{1,1}, pertrial_mpfr_noon{1,1});

        % Baseline firing rate statistics (Kolmogorov-Smirnov Test)
        [~, bsfr_ks_test_morn] = kstest(pertrial_bsfr_morn{1,1});
        [~, bsfr_ks_test_noon] = kstest(pertrial_bsfr_noon{1,1});
        % Baseline firing rate statistics (Unpaired T-Test)
        [~, bsfr_t_test] = ttest2(pertrial_bsfr_morn{1,1}, pertrial_bsfr_noon{1,1});
        % Baseline firing rate statistics (Wilcoxon rank sum test)
        [bsfr_wilcoxon, ~] = ranksum(pertrial_bsfr_morn{1,1}, pertrial_bsfr_noon{1,1});
        % Baseline firing rate percent change
        bsfr_perc = (avg_bsfr_noon - avg_bsfr_morn) / avg_bsfr_morn;
        % Baseline firing rate effect size (Cohen d)
        bsfr_cohen_d = Cohen_D(pertrial_bsfr_morn{1,1}, pertrial_bsfr_noon{1,1});

        % Depth of modulation statistics (Kolmogorov-Smirnov Test)
        [~, depth_ks_test_morn] = kstest(pertrial_depth_morn{1,1});
        [~, depth_ks_test_noon] = kstest(pertrial_depth_noon{1,1});
        % Depth of modulation statistics (Unpaired T-Test)
        [~, depth_t_test] = ttest2(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});
        % Depth of modulation statistics (Wilcoxon rank sum test)
        [depth_wilcoxon, ~] = ranksum(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});
        % Depth of modulation percent change
        depth_perc = (avg_depth_noon - avg_depth_morn) / avg_depth_morn;
        % Depth of modulation effect size (Cohen d)
        depth_cohen_d = Cohen_D(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});

        %% Add to the counter if above the depth minimum & SNR minimum
        cc = cc + 1;

        %% Put this info in the output arrays
        % Unit names
        all_unit_names{xx,1}(cc,1) = {unit_name};
        % Preferred Direction
        all_units_pref_dir{xx,1}(cc,1) = pref_dir;
        % Target Distance
        all_units_target{xx,1}(cc,1) = target_center;
        % Alignment Time
        alignment{xx,1}(cc,1) = alignment_time;
        % Post Spike Facilitation
        post_spike_facil{xx,1}(cc,1) = max_facil;

        % Morning baseline firing rate
        bsfr_morn{xx,1}(cc,1) = avg_bsfr_morn;
        % Afternoon baseline firing rate
        bsfr_noon{xx,1}(cc,1) = avg_bsfr_noon;
        % Morning depth
        depth_morn{xx,1}(cc,1) = avg_depth_morn;
        % Afternoon depth
        depth_noon{xx,1}(cc,1) = avg_depth_noon;

        % Morning baseline firing rate error
        bsfr_err_morn{xx,1}(cc,1) = err_bsfr_morn;
        % Morning depth error
        depth_err_morn{xx,1}(cc,1) = err_depth_morn;
        % Afternoon baseline firing rate error
        bsfr_err_noon{xx,1}(cc,1) = err_bsfr_noon;
        % Afternoon depth error
        depth_err_noon{xx,1}(cc,1) = err_depth_noon;

        % Modulation statistics
        mod_p_value_t_test_morn{xx,1}(cc,1) = mod_t_test_morn;
        mod_p_value_wilcoxon_morn{xx,1}(cc,1) = mod_wilcoxon_morn;
        mod_p_value_t_test_noon{xx,1}(cc,1) = mod_t_test_noon;
        mod_p_value_wilcoxon_noon{xx,1}(cc,1) = mod_wilcoxon_noon;

        % Baseline firing rate statistics
        bsfr_p_value_ks_test_morn{xx,1}(cc,1) = bsfr_ks_test_morn;
        bsfr_p_value_ks_test_noon{xx,1}(cc,1) = bsfr_ks_test_noon;
        bsfr_p_value_t_test{xx,1}(cc,1) = bsfr_t_test;
        bsfr_p_value_wilcoxon{xx,1}(cc,1) = bsfr_wilcoxon;
        bsfr_effect_perc{xx,1}(cc,1) = bsfr_perc;
        bsfr_effect_cohen_d{xx,1}(cc,1) = bsfr_cohen_d;

        % Depth of modulation statistics
        depth_p_value_ks_test_morn{xx,1}(cc,1) = depth_ks_test_morn;
        depth_p_value_ks_test_noon{xx,1}(cc,1) = depth_ks_test_noon;
        depth_p_value_t_test{xx,1}(cc,1) = depth_t_test;
        depth_p_value_wilcoxon{xx,1}(cc,1) = depth_wilcoxon;
        depth_effect_perc{xx,1}(cc,1) = depth_perc;
        depth_effect_cohen_d{xx,1}(cc,1) = depth_cohen_d;

        % Signal to noise ratio
        wave_sigtonoise{xx,1}(cc,1) = wave_sig_to_noise;
        nonlin_sigtonoise{xx,1}(cc,1) = nonlin_sig_to_noise;
        % Peak to peak amplitude
        wave_peaktopeak{xx,1}(cc,1) = wave_peak_to_peak;
        nonlin_peaktopeak{xx,1}(cc,1) = nonlin_peak_to_peak;
        % Wave shape t-test
        wave_p_value{xx,1}(cc,1) = wave_sort_metric;
        % Nonlinear energy t-test
        nonlin_p_value{xx,1}(cc,1) = nonlin_sort_metric;
        % Spike width
        spike_width{xx,1}(cc,1) = wave_spike_width;
        % Repolarization Time
        repol_time{xx,1}(cc,1) = wave_repol_time;
        % Fractional Contamination
        fract_contam{xx,1}(cc,1) = Fractional_Contam;

        % Number of targets in the preferred direction
        num_targets{xx,1}(cc,1) = pref_dir_targets;
  
    end % End the unit loop

    %% Create matrix for all units morning and afternoon

    % Create the excel matrix
    [morn_and_noon] = Multi_Session_Matrix(length(all_unit_names{xx,1}));

    % Assign the values to the excel matrix
    morn_and_noon.unit_names = char(all_unit_names{xx,1});
    morn_and_noon.bsfr_morn = bsfr_morn{xx,1};
    morn_and_noon.bsfr_noon = bsfr_noon{xx,1};
    morn_and_noon.depth_morn = depth_morn{xx,1};
    morn_and_noon.depth_noon = depth_noon{xx,1};
    morn_and_noon.bsfr_err_morn = bsfr_err_morn{xx,1};
    morn_and_noon.bsfr_err_noon = bsfr_err_noon{xx,1};
    morn_and_noon.depth_err_morn = depth_err_morn{xx,1};
    morn_and_noon.depth_err_noon = depth_err_noon{xx,1};
    morn_and_noon.bsfr_ks_test_morn = bsfr_p_value_ks_test_morn{xx,1};
    morn_and_noon.bsfr_ks_test_noon = bsfr_p_value_ks_test_noon{xx,1};
    morn_and_noon.bsfr_t_test = bsfr_p_value_t_test{xx,1};
    morn_and_noon.bsfr_wilcoxon = bsfr_p_value_wilcoxon{xx,1};
    morn_and_noon.bsfr_perc = bsfr_effect_perc{xx,1};
    morn_and_noon.bsfr_cohen_d = bsfr_effect_cohen_d{xx,1};
    morn_and_noon.mod_t_test_morn = mod_p_value_t_test_morn{xx,1};
    morn_and_noon.mod_wilcoxon_morn = mod_p_value_wilcoxon_morn{xx,1};
    morn_and_noon.mod_t_test_noon = mod_p_value_t_test_noon{xx,1};
    morn_and_noon.mod_wilcoxon_noon = mod_p_value_wilcoxon_noon{xx,1};
    morn_and_noon.depth_ks_test_morn = depth_p_value_ks_test_morn{xx,1};
    morn_and_noon.depth_ks_test_noon = depth_p_value_ks_test_noon{xx,1};
    morn_and_noon.depth_t_test = depth_p_value_t_test{xx,1};
    morn_and_noon.depth_wilcoxon = depth_p_value_wilcoxon{xx,1};
    morn_and_noon.depth_perc = depth_effect_perc{xx,1};
    morn_and_noon.depth_cohen_d = depth_effect_cohen_d{xx,1};
    morn_and_noon.pref_dir = all_units_pref_dir{xx,1};
    morn_and_noon.target = all_units_target{xx,1};
    morn_and_noon.alignment = alignment{xx,1};
    morn_and_noon.post_spike_facil = post_spike_facil{xx,1};
    morn_and_noon.wave_sigtonoise = wave_sigtonoise{xx,1};
    morn_and_noon.wave_peaktopeak = wave_peaktopeak{xx,1};
    morn_and_noon.wave_p_value = wave_p_value{xx,1};
    morn_and_noon.nonlin_sigtonoise = nonlin_sigtonoise{xx,1};
    morn_and_noon.nonlin_peaktopeak = nonlin_peaktopeak{xx,1};
    morn_and_noon.nonlin_p_value = nonlin_p_value{xx,1};
    morn_and_noon.spike_width = spike_width{xx,1};
    morn_and_noon.repol_time = repol_time{xx,1};
    morn_and_noon.fract_contam = fract_contam{xx,1};
    morn_and_noon.num_targets = num_targets{xx,1};
    morn_and_noon.rxn_time = NaN(height(morn_and_noon), 1);
    morn_and_noon.rxn_time(1) = rxn_time{xx,1};
    morn_and_noon.drug_dose_mg_per_kg = strings(height(morn_and_noon), 1);
    morn_and_noon.drug_dose_mg_per_kg(1) = char(Drug_Dose{xx,1});

    %% Save to Excel

    if isequal(Save_Excel, 1)

        % Define the file name
        filename = char(strcat(Dates{xx,1}, '_', Monkey, '_', ...
            Tasks{xx,1}, '_', Drug_Choice));

        % Save the file
        if ~exist(Save_Path, 'dir')
            mkdir(Save_Path);
        end
        writetable(morn_and_noon, strcat(Save_Path, filename, '.xlsx'))

    end

end % End the xds loop




