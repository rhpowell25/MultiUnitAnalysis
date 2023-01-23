function [bsfr_morn, bsfr_noon, depth_morn, depth_noon, bsfr_err_morn, mp_err_morn, bsfr_err_noon, mp_err_noon, ...
    depth_p_value_t_test, wave_p_value, nonlin_p_value, fract_contam, pref_dir, Drug_Dose, all_unit_names] = ...
    Multi_Session_Depth_Flex_Ext(Monkey, wrist_dir, Drug_Choice, Save_Excel)

%% Display the function being used
disp('Multi-Experiment Unit Results:');

%% Some of the analysis specifications

% What minimum depth of modulation? (#, or NaN)
Depth_Minimum = NaN;

% What minimum signal to noise ration? (#, or NaN)
SNR_Minimum = 5;

% Which targets do you want the mnovement phase firing rate calculated from? ('Max', 'Min', 'All')
tgt_mpfr = 'Max';

Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel Data\', wrist_dir, '\', tgt_mpfr, '_Targets\');

% Load the file information
[Dates, Tasks, Drug_Dose] = File_Details(Monkey, Drug_Choice);

%% Build the output arrays

% Morning baseline firing rate
bsfr_morn = struct([]);
% Morning depth
depth_morn = struct([]);
% Afternoon baseline firing rate
bsfr_noon = struct([]);
% Afternoon depth
depth_noon = struct([]);
    
% Morning error
bsfr_err_morn = struct([]);
mp_err_morn = struct([]);
% Afternoon error
bsfr_err_noon = struct([]);
mp_err_noon = struct([]);

% Depth of modulation t-test
depth_p_value_t_test = struct([]);
depth_p_value_wilcoxon = struct([]);

% Wave shape t-test
wave_p_value = struct([]);
% Nonlinear energy t-test
nonlin_p_value = struct([]);
% Fractional Contamination
fract_contam = struct([]);

% Preferred Direction
pref_dir = struct([]);

% Unit Names
all_unit_names = struct([]);

% Number of targets
num_targets = struct([]);

%% Loop through the different experiments
for xx = 1:length(Dates)
    
    % Load the relevant xds file
    xds_morn = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, 'Noon');

    % Process the xds files
    Match_The_Targets = 0;

    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);

    % Define the target direction
    if strcmp(wrist_dir, 'Flex')
        if strcmp(xds_morn.meta.hand, 'Left')
            selec_dir = 0;
        end
        if strcmp(xds_morn.meta.hand, 'Right')
            selec_dir = 180;
        end
    end
    if strcmp(wrist_dir, 'Ext')
        if strcmp(xds_morn.meta.hand, 'Left')
            selec_dir = 180;
        end
        if strcmp(xds_morn.meta.hand, 'Right')
            selec_dir = 0;
        end
    end

    % Extract the target directions & centers
    [target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
    [target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

    % Check to see if both sessions use a consistent number of targets
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

    % Skip the file if no trials in the selected direction
    if ~any(ismember(target_dirs_morn, selec_dir))
        continue
    end

    % Set the succesful unit counter
    cc = 0;

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

        %% Get the baseline firing rates
        [~, ~, ~, pertrial_bsfr_morn] = ... 
            BaselineFiringRate(xds_morn, unit_name);
        [~, ~, ~, pertrial_bsfr_noon] = ... 
            BaselineFiringRate(xds_noon, unit_name);

        %% Get the movement phase firing rates

        event = 'window_trial_gocue';

        disp('Morning')
        [~, ~, ~, pertrial_mpfr_morn] = ...
            EventFiringRate(xds_morn, unit_name, event, tgt_mpfr);
        disp('Afternoon')
        [~, ~, ~, pertrial_mpfr_noon] = ... 
            EventFiringRate(xds_noon, unit_name, event, tgt_mpfr);

        %% Only use the info of target centers conserved between morn & noon
        if ~all(Matching_Idxs_Morn) || ~all(Matching_Idxs_Noon)
            disp('Uneven Targets Between Morning & Afternoon');
            pertrial_mpfr_morn = pertrial_mpfr_morn(Matching_Idxs_Morn);
            pertrial_mpfr_noon = pertrial_mpfr_noon(Matching_Idxs_Noon);
        end

        %% Calculate the depth of modulation per trial
        pertrial_depth_morn = struct([]);
        pertrial_depth_noon = struct([]);
        for ii = 1:length(pertrial_mpfr_morn)
            pertrial_depth_morn{ii,1} = pertrial_mpfr_morn{ii,1} - mean(pertrial_bsfr_morn{1,1});
            pertrial_depth_noon{ii,1} = pertrial_mpfr_noon{ii,1} - mean(pertrial_bsfr_noon{1,1});
        end
        
        %% Finding the average & standard error of the depth of modulation & baseline firing rates

        avg_bsfr_morn = zeros(length(pertrial_mpfr_morn),1);
        avg_depth_morn = zeros(length(pertrial_mpfr_morn),1);
        err_bsfr_morn = zeros(length(pertrial_mpfr_morn),1);
        err_depth_morn = zeros(length(pertrial_mpfr_morn),1);
        avg_bsfr_noon = zeros(length(pertrial_mpfr_noon),1);
        avg_depth_noon = zeros(length(pertrial_mpfr_noon),1);
        err_depth_noon = zeros(length(pertrial_mpfr_noon),1);
        err_bsfr_noon = zeros(length(pertrial_mpfr_noon),1);
        for ii = 1:length(avg_depth_morn)
            avg_bsfr_morn(ii,1) = mean(pertrial_bsfr_morn{1,1});
            avg_depth_morn(ii,1) = mean(pertrial_depth_morn{ii,1});
            err_bsfr_morn(ii,1) = std(pertrial_bsfr_morn{1,1}) / ...
                sqrt(length(pertrial_bsfr_morn{1,1}));
            err_depth_morn(ii,1) = std(pertrial_depth_morn{ii,1}) / ...
                sqrt(length(pertrial_depth_morn{ii,1}));
            avg_bsfr_noon(ii,1) = mean(pertrial_bsfr_noon{1,1});
            avg_depth_noon(ii,1) = mean(pertrial_depth_noon{ii,1});
            err_bsfr_noon(ii,1) = std(pertrial_bsfr_noon{1,1}) / ...
                sqrt(length(pertrial_bsfr_noon{1,1}));
            err_depth_noon(ii,1) = std(pertrial_depth_noon{ii,1}) / ...
                sqrt(length(pertrial_depth_noon{ii,1}));
        end

        %% Only look at the selected direction

        selec_dir_max_tgt_morn = max(target_centers_morn(target_dirs_morn == selec_dir));
        selec_dir_min_tgt_morn = min(target_centers_morn(target_dirs_morn == selec_dir));
        if isempty(selec_dir_max_tgt_morn) || isempty(selec_dir_min_tgt_morn)
            disp('No targets in the preferred direction!')
            continue
        end

        selec_dir_max_tgt_noon = max(target_centers_noon(target_dirs_noon == selec_dir));
        selec_dir_min_tgt_noon = min(target_centers_noon(target_dirs_noon == selec_dir));
        if isempty(selec_dir_max_tgt_noon) || isempty(selec_dir_min_tgt_noon)
            disp('No targets in the preferred direction!')
            continue
        end
        
        %% Look at the maximum or minimum targets if not using all targets
        if strcmp(tgt_mpfr, 'Max')
            tgt_bsfr_morn = avg_bsfr_morn(target_centers_morn == selec_dir_max_tgt_morn);
            tgt_bsfr_noon = avg_bsfr_noon(target_centers_noon == selec_dir_max_tgt_noon);
    
            tgt_depth_morn = avg_depth_morn(target_centers_morn == selec_dir_max_tgt_morn);
            tgt_depth_noon = avg_depth_noon(target_centers_noon == selec_dir_max_tgt_noon);
    
            tgt_pertrial_depth_morn = pertrial_depth_morn(target_centers_morn == selec_dir_max_tgt_morn);
            tgt_pertrial_depth_noon = pertrial_depth_noon(target_centers_noon == selec_dir_max_tgt_noon);
    
            tgt_err_bsfr_morn = err_bsfr_morn(target_centers_morn == selec_dir_max_tgt_morn);
            tgt_err_bsfr_noon = err_bsfr_noon(target_centers_noon == selec_dir_max_tgt_noon);
    
            tgt_err_depth_morn = err_depth_morn(target_centers_morn == selec_dir_max_tgt_morn);
            tgt_err_depth_noon = err_depth_noon(target_centers_noon == selec_dir_max_tgt_noon);
    
            tgt_dir_morn = target_dirs_morn(target_centers_morn == selec_dir_max_tgt_morn);
            tgt_dir_noon = target_dirs_noon(target_centers_noon == selec_dir_max_tgt_noon);
        end

        if strcmp(tgt_mpfr, 'Min')
            tgt_bsfr_morn = avg_bsfr_morn(target_centers_morn == selec_dir_min_tgt_morn);
            tgt_bsfr_noon = avg_bsfr_noon(target_centers_noon == selec_dir_min_tgt_noon);
    
            tgt_depth_morn = avg_depth_morn(target_centers_morn == selec_dir_min_tgt_morn);
            tgt_depth_noon = avg_depth_noon(target_centers_noon == selec_dir_min_tgt_noon);
    
            tgt_pertrial_depth_morn = pertrial_depth_morn(target_centers_morn == selec_dir_min_tgt_morn);
            tgt_pertrial_depth_noon = pertrial_depth_noon(target_centers_noon == selec_dir_min_tgt_noon);
    
            tgt_err_bsfr_morn = err_bsfr_morn(target_centers_morn == selec_dir_min_tgt_morn);
            tgt_err_bsfr_noon = err_bsfr_noon(target_centers_noon == selec_dir_min_tgt_noon);
    
            tgt_err_depth_morn = err_depth_morn(target_centers_morn == selec_dir_min_tgt_morn);
            tgt_err_depth_noon = err_depth_noon(target_centers_noon == selec_dir_min_tgt_noon);
    
            tgt_dir_morn = target_dirs_morn(target_centers_morn == selec_dir_min_tgt_morn);
            tgt_dir_noon = target_dirs_noon(target_centers_noon == selec_dir_min_tgt_noon);
        end

        %% Find the number of targets in the selected direction
        num_targets_morn = length(target_centers_morn(target_dirs_morn == selec_dir));
        num_targets_noon = length(target_centers_morn(target_dirs_morn == selec_dir));
        % Confirm the # of targets in the morning & afternoon are equal
        if isequal(num_targets_morn, num_targets_noon)
            pref_dir_targets = num_targets_morn;
        else
            pref_dir_targets = NaN;
            disp('Unequal number of targets!')
        end
  
        %% Baseline firing rate & depth of modulation in the selected direction
        avg_bsfr_morn = tgt_bsfr_morn(tgt_dir_morn == selec_dir);
        avg_bsfr_noon = tgt_bsfr_noon(tgt_dir_noon == selec_dir);

        pertrial_depth_morn = tgt_pertrial_depth_morn(tgt_dir_morn == selec_dir);
        pertrial_depth_noon = tgt_pertrial_depth_noon(tgt_dir_noon == selec_dir);

        avg_depth_morn = tgt_depth_morn(tgt_dir_morn == selec_dir);
        avg_depth_noon = tgt_depth_noon(tgt_dir_noon == selec_dir);

        %% Skip the unit if the depths of modulation is below the defined minimum
        if ~isnan(Depth_Minimum)
            if avg_depth_morn <= Depth_Minimum && avg_depth_noon <= Depth_Minimum
                continue
            end
        end

        %% Skip the unit if the SNR is below the defined minimum

        if ~isnan(SNR_Minimum)
            [SNR_matrix_morn] = SignalToNoise(xds_morn, unit_name, 0);
            SNR_morn = SNR_matrix_morn{2,1};
            [SNR_matrix_noon] = SignalToNoise(xds_noon, unit_name, 0);
            SNR_noon = SNR_matrix_noon{2,1};
            if SNR_morn <= SNR_Minimum || SNR_noon <= SNR_Minimum
                continue
            end
        end

        %% Standard error of the baseline firing rate & depth of modulation in the selected direction
        err_bsfr_morn = tgt_err_bsfr_morn(tgt_dir_morn == selec_dir);
        err_bsfr_noon = tgt_err_bsfr_noon(tgt_dir_noon == selec_dir);

        err_depth_morn = tgt_err_depth_morn(tgt_dir_morn == selec_dir);
        err_depth_noon = tgt_err_depth_noon(tgt_dir_noon == selec_dir);

        %% Check if the unit's depth of modulation changed significantly

        % Depth of modulation statistics (Unpaired T-Test)
        [~, depth_t_test] = ttest2(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});

        % Depth of modulation statistics (Wilcoxon rank sum test)
        [depth_wilcoxon, ~] = ranksum(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});

        %% Check if the unit is well sorted
        

        [~, wave_sort_metric] = WaveShapes_Morn_v_Noon(xds_morn, xds_noon, unit_name, 0, 0);
        if isnan(wave_sort_metric)
            % If the unit doesn't exist or has less than 1000 spikes
            continue
        end
        [~, nonlin_sort_metric] = NonLinearEnergy_Morn_v_Noon(xds_morn, xds_noon, unit_name, 0, 0);


        [Fractional_Contam, ~] = InterstimulusInterval_Morn_v_Noon(xds_morn, xds_noon, unit_name, 0, 0);

        [fract_contam_morn, ~] = InterstimulusInterval(xds_morn, unit_name, 0, 0);
        [fract_contam_noon, ~] = InterstimulusInterval(xds_noon, unit_name, 0, 0);

        % Skip the unit if the fract_contam changes between sessions
        if abs(fract_contam_noon - fract_contam_morn) >= 0.5
            continue
        end

        %% Add to the counter if above the depth minimum & SNR minimum
        cc = cc + 1;

        %% Put this info in the output arrays
        % Unit names
        all_unit_names{xx,1}(cc,1) = {unit_name};
        % Preferred Direction
        pref_dir{xx,1}(cc,1) = selec_dir;

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
        mp_err_morn{xx,1}(cc,1) = err_depth_morn;
        % Afternoon baseline firing rate error
        bsfr_err_noon{xx,1}(cc,1) = err_bsfr_noon;
        % Afternoon depth error
        mp_err_noon{xx,1}(cc,1) = err_depth_noon;

        % Depth of modulation statistics
        depth_p_value_t_test{xx,1}(cc,1) = depth_t_test;
        depth_p_value_wilcoxon{xx,1}(cc,1) = depth_wilcoxon;

        % Wave shape t-test
        wave_p_value{xx,1}(cc,1) = wave_sort_metric;
        % Nonlinear energy t-test
        nonlin_p_value{xx,1}(cc,1) = nonlin_sort_metric;
        % Fractional Contamination
        fract_contam{xx,1}(cc,1) = Fractional_Contam;

        % Number of targets in the preferred direction
        num_targets{xx,1}(cc,1) = pref_dir_targets;
  
    end % End the unit loop

    %% Create matrix for all units morning and afternoon

    % Create the excel matrix
    morn_and_noon = cell(length(all_unit_names{xx,1}) + 1, 17);
    % Define the excel matrix headers
    morn_and_noon{1,1} = 'unit_names';
    morn_and_noon{1,2} = 'pref_dir';
    morn_and_noon{1,3} = 'bsfr_morn';
    morn_and_noon{1,4} = 'bsfr_err_morn';
    morn_and_noon{1,5} = 'mp_err_morn';
    morn_and_noon{1,6} = 'depth_morn';
    morn_and_noon{1,7} = 'bsfr_noon';
    morn_and_noon{1,8} = 'bsfr_err_noon';
    morn_and_noon{1,9} = 'mp_err_noon';
    morn_and_noon{1,10} = 'depth_noon';
    morn_and_noon{1,11} = 'depth_t_test';
    morn_and_noon{1,12} = 'depth_wilcoxon';
    morn_and_noon{1,13} = 'wave_p_value';
    morn_and_noon{1,14} = 'nonlin_p_value';
    morn_and_noon{1,15} = 'fract_contam';
    morn_and_noon{1,16} = 'num_targets';
    morn_and_noon{1,17} = 'drug_dose_mg_per_kg';

    % Assign the values to the excel matrix
    uu = 2;
    morn_and_noon{2,17} = Drug_Dose{xx,1};
    for N = 1:length(all_unit_names{xx,1})
        morn_and_noon{uu,1} = char(all_unit_names{xx,1}(N,1));
        morn_and_noon{uu,2} = pref_dir{xx,1}(N,1);
        morn_and_noon{uu,3} = bsfr_morn{xx,1}(N,1);
        morn_and_noon{uu,4} = bsfr_err_morn{xx,1}(N,1);
        morn_and_noon{uu,5} = mp_err_morn{xx,1}(N,1);
        morn_and_noon{uu,6} = depth_morn{xx,1}(N,1);
        morn_and_noon{uu,7} = bsfr_noon{xx,1}(N,1);
        morn_and_noon{uu,8} = bsfr_err_noon{xx,1}(N,1);
        morn_and_noon{uu,9} = mp_err_noon{xx,1}(N,1);
        morn_and_noon{uu,10} = depth_noon{xx,1}(N,1);
        morn_and_noon{uu,11} = depth_p_value_t_test{xx,1}(N,1);
        morn_and_noon{uu,12} = depth_p_value_wilcoxon{xx,1}(N,1);
        morn_and_noon{uu,13} = wave_p_value{xx,1}(N,1);
        morn_and_noon{uu,14} = nonlin_p_value{xx,1}(N,1);
        morn_and_noon{uu,15} = fract_contam{xx,1}(N,1);
        morn_and_noon{uu,16} = num_targets{xx,1}(N,1);
        uu = uu + 1;
    end

    %% Save to Excel

    if isequal(Save_Excel, 1)

        % Define the file name
        Monkey_Name = xds_morn.meta.monkey;
        filename = char(strcat(Dates{xx,1}, '_', Monkey_Name, '_', ...
            Tasks{xx,1}, '_', Drug_Choice));

        % Save the file
        if ~exist(Save_Path, 'dir')
            mkdir(Save_Path);
        end
        writecell(morn_and_noon, strcat(Save_Path, filename, '.xlsx'))

    end

end % End the xds loop




