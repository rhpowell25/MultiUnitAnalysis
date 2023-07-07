function Multi_Session_Split(Monkey, event, Drug_Choice, Save_Excel)

%% Display the function being used
disp('Multi-Experiment Unit Results:');

%% Some of the analysis specifications

% How do you want to split the afternoon triasl ('Random', 'Half')
Split_Method = 'Random';

% Which targets do you want the mnovement phase firing rate calculated from? ('Max', 'Min', 'All')
tgt_mpfr = 'Max';

% Load the sorted files
Sorted = 1;

if isequal(Sorted, 1)
    Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '_1\Sorted\');
else
    Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '_1\Unsorted\');
end

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug_Choice);

%% Build the output arrays

% Afternoon depth
depth_noon = struct([]);

% Afternoon ramp phase firing rate
ramp_noon = struct([]);

% Afternoon TgtHold phase firing rate
TgtHold_noon = struct([]);

% Afternoon error
ramp_err_noon = struct([]);
TgtHold_err_noon = struct([]);
mp_err_noon = struct([]);

% Modulation statistics
mod_p_value_t_test_noon = struct([]);
mod_p_value_wilcoxon_noon = struct([]);

% Ramp phase statistics
ramp_p_value_t_test = struct([]);
ramp_p_value_wilcoxon = struct([]);
ramp_effect_perc = struct([]);
ramp_effect_cohen_d = struct([]);
% TgtHold  phase statistics
TgtHold_p_value_t_test = struct([]);
TgtHold_p_value_wilcoxon = struct([]);
TgtHold_effect_perc = struct([]);
TgtHold_effect_cohen_d = struct([]);
% Depth of modulation statistics
depth_p_value_t_test = struct([]);
depth_p_value_wilcoxon = struct([]);
depth_effect_perc = struct([]);
depth_effect_cohen_d = struct([]);

%% Loop through the different experiments
for xx = 1:length(Dates)
    
    % Load the relevant xds file
    xds_morn = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Noon');

    % Process the xds files
    Match_The_Targets = 0;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);

    % Load the excel file
    [xds_excel] = Load_Excel(Monkey, Dates{xx,1}, Tasks{xx,1});

    % Find the number of succesful trials
    rewarded_idx = Find_Max_Indexes(xds_noon);
    
    % Divide the number of trials by two
    split_length = round(length(rewarded_idx) / 2);
    
    if strcmp(Split_Method, 'Random')
        first_split_idxs = randperm(length(rewarded_idx), split_length);
    elseif strcmp(Split_Method, 'Half')
        first_split_idxs = 1:split_length;
    end

    %% Loop twice for each of the trial splits
    for pp = 1:2

        if isequal(pp, 1)
            split_idxs = first_split_idxs;
        else
            if strcmp(Split_Method, 'Random')
                temp = 1:length(rewarded_idx);
                split_idxs = setdiff(temp, first_split_idxs);
            elseif strcmp(Split_Method, 'Half')
                split_idxs = split_length + 1:length(rewarded_idx);
            end
        end

        %% Loop through the units  
        for jj = 1:length(xds_excel.unit_names)
    
            %% Define the unit
            unit_name = char(xds_excel.unit_names(jj));
    
            %% Get the baseline firing rates
            [avg_bsfr_morn, ~, ~, ~] = ... 
                BaselineFiringRate(xds_morn, unit_name);
            [avg_bsfr_noon, ~, ~, pertrial_bsfr_noon] = ... 
                BaselineFiringRate(xds_noon, unit_name);
    
            %% Get the ramp phase firing rates
            [avg_ramp_morn, ~, ~, pertrial_ramp_morn] = ... 
                RampFiringRate(xds_morn, unit_name);
            [~, ~, ~, pertrial_ramp_noon] = ... 
                RampFiringRate(xds_noon, unit_name);
    
            %% Get the target hold firing rates
            [avg_TgtHold_morn, ~, ~, pertrial_TgtHold_morn] = ... 
                TgtHoldFiringRate(xds_morn, unit_name);
            [~, ~, ~, pertrial_TgtHold_noon] = ... 
                TgtHoldFiringRate(xds_noon, unit_name);
    
            %% Get the movement phase firing rates
    
            [~, ~, pertrial_mpfr_morn] = ...
                EventPeakFiringRate(xds_morn, unit_name, event);
            [~, ~, pertrial_mpfr_noon] = ... 
                EventPeakFiringRate(xds_noon, unit_name, event);
    
            % Extract the target directions & centers
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
                pertrial_mpfr_morn = pertrial_mpfr_morn(Matching_Idxs_Morn);
                pertrial_mpfr_noon = pertrial_mpfr_noon(Matching_Idxs_Noon);
            end
    
            %% Calculate the depth of modulation per trial
            pertrial_depth_morn = struct([]);
            pertrial_depth_noon = struct([]);
            for ii = 1:length(pertrial_mpfr_morn)
                pertrial_depth_morn{ii,1} = pertrial_mpfr_morn{ii,1} - avg_bsfr_morn(1,1);
                pertrial_depth_noon{ii,1} = pertrial_mpfr_noon{ii,1} - avg_bsfr_noon(1,1);
            end
            
            %% Find the mean & standard error of the depth of modulation
    
            avg_depth_morn = zeros(length(pertrial_mpfr_morn),1);
            err_depth_morn = zeros(length(pertrial_mpfr_morn),1);
            for ii = 1:length(avg_depth_morn)
                avg_depth_morn(ii,1) = mean(pertrial_depth_morn{ii,1});
                err_depth_morn(ii,1) = std(pertrial_depth_morn{ii,1}) / ...
                    sqrt(length(pertrial_depth_morn{ii,1}));
            end
    
            %% Only look at the preferred direction
            [pref_dir] = PreferredDirection_Morn_v_Noon(xds_morn, xds_noon, unit_name, event, tgt_mpfr);
    
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
    
            %% Look at the maximum or minimum targets if not using all targets
    
            tgt_ramp_morn = avg_ramp_morn(target_centers_morn == pref_dir_tgt_morn);
            tgt_TgtHold_morn = avg_TgtHold_morn(target_centers_morn == pref_dir_tgt_morn);
            tgt_depth_morn = avg_depth_morn(target_centers_morn == pref_dir_tgt_morn);
    
            tgt_pertrial_ramp_morn = pertrial_ramp_morn(target_centers_morn == pref_dir_tgt_morn);
            tgt_pertrial_ramp_noon = pertrial_ramp_noon(target_centers_noon == pref_dir_tgt_noon);
    
            tgt_pertrial_TgtHold_morn = pertrial_TgtHold_morn(target_centers_morn == pref_dir_tgt_morn);
            tgt_pertrial_TgtHold_noon = pertrial_TgtHold_noon(target_centers_noon == pref_dir_tgt_noon);
    
            tgt_pertrial_mpfr_noon = pertrial_mpfr_noon(target_centers_noon == pref_dir_tgt_noon);
    
            tgt_pertrial_depth_morn = pertrial_depth_morn(target_centers_morn == pref_dir_tgt_morn);
            tgt_pertrial_depth_noon = pertrial_depth_noon(target_centers_noon == pref_dir_tgt_noon);
    
            tgt_dir_morn = target_dirs_morn(target_centers_morn == pref_dir_tgt_morn);
            tgt_dir_noon = target_dirs_noon(target_centers_noon == pref_dir_tgt_noon);
      
            %% Phasic firing rates & depth of modulation in the preferred direction
    
            avg_ramp_morn = tgt_ramp_morn(tgt_dir_morn == pref_dir);
    
            avg_TgtHold_morn = tgt_TgtHold_morn(tgt_dir_morn == pref_dir);
    
            pertrial_ramp_morn = tgt_pertrial_ramp_morn(tgt_dir_morn == pref_dir);
            pertrial_ramp_noon = tgt_pertrial_ramp_noon(tgt_dir_noon == pref_dir);
    
            pertrial_TgtHold_morn = tgt_pertrial_TgtHold_morn(tgt_dir_morn == pref_dir);
            pertrial_TgtHold_noon = tgt_pertrial_TgtHold_noon(tgt_dir_noon == pref_dir);
    
            pertrial_mpfr_noon = tgt_pertrial_mpfr_noon(tgt_dir_noon == pref_dir);
    
            pertrial_depth_morn = tgt_pertrial_depth_morn(tgt_dir_morn == pref_dir);
            pertrial_depth_noon = tgt_pertrial_depth_noon(tgt_dir_noon == pref_dir);
    
            avg_depth_morn = tgt_depth_morn(tgt_dir_morn == pref_dir);

            %% Subsample the trials according to the trial split

            pertrial_ramp_noon{1,1} = pertrial_ramp_noon{1,1}(split_idxs);
            avg_ramp_noon = mean(pertrial_ramp_noon{1,1});
            err_ramp_noon = std(pertrial_ramp_noon{1,1}) / ...
                    sqrt(length(pertrial_ramp_noon{1,1}));
            pertrial_TgtHold_noon{1,1} = pertrial_TgtHold_noon{1,1}(split_idxs);
            avg_TgtHold_noon = mean(pertrial_TgtHold_noon{1,1});
            err_TgtHold_noon = std(pertrial_TgtHold_noon{1,1}) / ...
                    sqrt(length(pertrial_TgtHold_noon{1,1}));
            pertrial_mpfr_noon{1,1} = pertrial_mpfr_noon{1,1}(split_idxs);
            pertrial_depth_noon{1,1} = pertrial_depth_noon{1,1}(split_idxs);
            avg_depth_noon = mean(pertrial_depth_noon{1,1});
            err_depth_noon = std(pertrial_depth_noon{1,1}) / ...
                    sqrt(length(pertrial_depth_noon{1,1}));
    
            %% Check if the unit's depth of modulation changed significantly
    
            % Modulation statistics & effect sizes
            [~, mod_t_test_noon] = ttest2(pertrial_bsfr_noon{1,1}, pertrial_mpfr_noon{1,1});
            [mod_wilcoxon_noon, ~] = ranksum(pertrial_bsfr_noon{1,1}, pertrial_mpfr_noon{1,1});
    
            % Ramp phase statistics (Unpaired T-Test)
            [~, ramp_t_test] = ttest2(pertrial_ramp_morn{1,1}, pertrial_ramp_noon{1,1});
            % Ramp phase statistics (Wilcoxon rank sum test)
            try
                [ramp_wilcoxon, ~] = ranksum(pertrial_ramp_morn{1,1}, pertrial_ramp_noon{1,1});
            catch
                ramp_wilcoxon = NaN;
            end
            % Ramp phase percent change
            ramp_perc = (avg_ramp_noon - avg_ramp_morn) / avg_ramp_morn;
            % Ramp phase effect size (Cohen d)
            ramp_cohen_d = Cohen_D(pertrial_ramp_morn{1,1}, pertrial_ramp_noon{1,1});
    
            % TgtHold phase statistics (Unpaired T-Test)
            [~, TgtHold_t_test] = ttest2(pertrial_TgtHold_morn{1,1}, pertrial_TgtHold_noon{1,1});
            % TgtHold phase statistics (Wilcoxon rank sum test)
            [TgtHold_wilcoxon, ~] = ranksum(pertrial_TgtHold_morn{1,1}, pertrial_TgtHold_noon{1,1});
            % TgtHold phase percent change
            TgtHold_perc = (avg_TgtHold_noon - avg_TgtHold_morn) / avg_TgtHold_morn;
            % TgtHold phase effect size (Cohen d)
            TgtHold_cohen_d = Cohen_D(pertrial_TgtHold_morn{1,1}, pertrial_TgtHold_noon{1,1});
    
            % Depth of modulation statistics (Unpaired T-Test)
            [~, depth_t_test] = ttest2(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});
            % Depth of modulation statistics (Wilcoxon rank sum test)
            [depth_wilcoxon, ~] = ranksum(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});
            % Depth of modulation percent change
            depth_perc = (avg_depth_noon - avg_depth_morn) / avg_depth_morn;
            % Depth of modulation effect size (Cohen d)
            depth_cohen_d = Cohen_D(pertrial_depth_morn{1,1}, pertrial_depth_noon{1,1});
    
            %% Put this info in the output arrays
    
            % Afternoon ramp firing rate
            ramp_noon{xx,1}(jj,1) = avg_ramp_noon;
    
            % Afternoon TgtHold firing rate
            TgtHold_noon{xx,1}(jj,1) = avg_TgtHold_noon;
            % Afternoon depth
            depth_noon{xx,1}(jj,1) = avg_depth_noon;
    
            % Afternoon ramp firing rate error
            ramp_err_noon{xx,1}(jj,1) = err_ramp_noon;
            % Afternoon TgtHold firing rate error
            TgtHold_err_noon{xx,1}(jj,1) = err_TgtHold_noon;
            % Afternoon depth error
            mp_err_noon{xx,1}(jj,1) = err_depth_noon;
    
            % Modulation statistics
            mod_p_value_t_test_noon{xx,1}(jj,1) = mod_t_test_noon;
            mod_p_value_wilcoxon_noon{xx,1}(jj,1) = mod_wilcoxon_noon;
    
            % Ramp phase statistics
            ramp_p_value_t_test{xx,1}(jj,1) = ramp_t_test;
            ramp_p_value_wilcoxon{xx,1}(jj,1) = ramp_wilcoxon;
            ramp_effect_perc{xx,1}(jj,1) = ramp_perc;
            ramp_effect_cohen_d{xx,1}(jj,1) = ramp_cohen_d;
    
            % TgtHold phase statistics
            TgtHold_p_value_t_test{xx,1}(jj,1) = TgtHold_t_test;
            TgtHold_p_value_wilcoxon{xx,1}(jj,1) = TgtHold_wilcoxon;
            TgtHold_effect_perc{xx,1}(jj,1) = TgtHold_perc;
            TgtHold_effect_cohen_d{xx,1}(jj,1) = TgtHold_cohen_d;
    
            % Depth of modulation statistics
            depth_p_value_t_test{xx,1}(jj,1) = depth_t_test;
            depth_p_value_wilcoxon{xx,1}(jj,1) = depth_wilcoxon;
            depth_effect_perc{xx,1}(jj,1) = depth_perc;
            depth_effect_cohen_d{xx,1}(jj,1) = depth_cohen_d;
      
        end % End the unit loop
    
        %% Replace the matrix for all units morning and afternoon
    
        % Assign the values to the excel matrix
        xds_excel.depth_noon = depth_noon{xx,1};
        xds_excel.mp_err_noon = mp_err_noon{xx,1};
        xds_excel.mod_t_test_noon = mod_p_value_t_test_noon{xx,1};
        xds_excel.mod_wilcoxon_noon = mod_p_value_wilcoxon_noon{xx,1};
        xds_excel.depth_t_test = depth_p_value_t_test{xx,1};
        xds_excel.depth_wilcoxon = depth_p_value_wilcoxon{xx,1};
        xds_excel.depth_perc = depth_effect_perc{xx,1};
        xds_excel.depth_cohen_d = depth_effect_cohen_d{xx,1};
        xds_excel.ramp_noon = ramp_noon{xx,1};
        xds_excel.ramp_err_noon = ramp_err_noon{xx,1};
        xds_excel.TgtHold_noon = TgtHold_noon{xx,1};
        xds_excel.TgtHold_err_noon = TgtHold_err_noon{xx,1};
        xds_excel.ramp_t_test = ramp_p_value_t_test{xx,1};
        xds_excel.ramp_wilcoxon = ramp_p_value_wilcoxon{xx,1};
        xds_excel.ramp_perc = ramp_effect_perc{xx,1};
        xds_excel.ramp_cohen_d = ramp_effect_cohen_d{xx,1};
        xds_excel.TgtHold_t_test = TgtHold_p_value_t_test{xx,1};
        xds_excel.TgtHold_wilcoxon = TgtHold_p_value_wilcoxon{xx,1};
        xds_excel.TgtHold_perc = TgtHold_effect_perc{xx,1};
        xds_excel.TgtHold_cohen_d = TgtHold_effect_cohen_d{xx,1};
    
        %% Save to Excel
    
        if isequal(Save_Excel, 1)
    
            % Define the file name
            filename = char(strcat(Dates{xx,1}, '_', Monkey, '_', ...
                Tasks{xx,1}, '_', Drug_Choice));
    
            % Save the file
            if ~exist(Save_Path, 'dir')
                mkdir(Save_Path);
            end
            writetable(xds_excel, strcat(Save_Path, filename, '_', ...
                Split_Method, '_', num2str(pp), '.xlsx'))
    
        end
    
    end % End of the split loop

end % End the xds loop




