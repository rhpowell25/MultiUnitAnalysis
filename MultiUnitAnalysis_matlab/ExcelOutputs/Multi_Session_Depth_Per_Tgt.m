function [depth, pref_dir, all_unit_names] = ...
    Multi_Session_Depth_Per_Tgt(Monkey, Drug_Choice, Save_Excel)

%% Display the function being used
disp('Multi-Experiment Unit Results:');

%% Some of the analysis specifications

% What minimum signal to noise ration? (#, or NaN)
SNR_Minimum = 5;

Save_Path = 'C:\Users\rhpow\Documents\Work\Northwestern\Excel Data\Noon_';

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug_Choice);

%% Build the output arrays

% Unit Names
all_unit_names = struct([]);

% Depth of modulation
depth = struct([]);

% Target center
Tgt_center = struct([]);

% Preferred Direction
units_pref_dir = struct([]);

%% Loop through the different experiments
for xx = 1:length(Dates)
    
    % Load the relevant xds file
    [~, xds, ~] = Load_XDS(Dates{xx,1}, Tasks{xx,1}, 0);

    % Extract the target directions & centers
    [target_dirs, target_centers] = Identify_Targets(xds);

    % Set the succesful unit counter
    cc = 1;

    %% Loop through the units  
    for uu = 1:length(xds.unit_names)

        %% Skip the unit if the SNR is below the defined minimum
        if ~isnan(SNR_Minimum)
            [SNR_matrix] = SignalToNoise(xds, char(xds.unit_names(uu)), 0);
            SNR = SNR_matrix{2,1};
            if SNR <= SNR_Minimum
                continue
            end
        end

        %% Get the baseline firing rates
        [~, ~, ~, ~, ~, pertrial_bsfr] = ... 
            BaselineFiringRate(xds, char(xds.unit_names(uu)));

        %% Loop through the target directions
        for dd = 1:length(target_dirs)

            % Find the number of targets in that direction
            num_targets = length(target_centers(target_dirs == target_dirs(dd)));

            % Define the mpfr based on the target
            if isequal(num_targets, 1)
                tgt_mpfr = 'Max';
            end

            if num_targets >= 2
                tgts_per_dir = target_centers(target_dirs == target_dirs(dd));
                if isequal(target_centers(dd), max(tgts_per_dir))
                    tgt_mpfr = 'Max';
                end
                if isequal(target_centers(dd), min(tgts_per_dir))
                    tgt_mpfr = 'Min';
                end
                % If there are three targets
                if rem(length(tgts_per_dir), 2) == 1
                    if isequal(target_centers(dd), tgts_per_dir(ceil(numel(tgts_per_dir)/2)))
                        tgt_mpfr = 'Mid';
                    end
                end
            end

            %% Only look at the preferred direction
            pref_dir = GoCuePreferredDirection(xds, char(xds.unit_names(uu)), tgt_mpfr);

            %% Continue if the preferred direction is mismatched
            if ~isequal(pref_dir, target_dirs(dd))
                continue
            end

            %% Get the movement phase firing rates
            [~, ~, ~, ~, ~, pertrial_mpfr] = ... 
                WindowTrialGoCueFiringRate(xds, char(xds.unit_names(uu)), tgt_mpfr);

            %% Depth of modulation in the preferred direction
            pertrial_mpfr = pertrial_mpfr(target_dirs == pref_dir);
            
            %% Calculate the depth of modulation per trial
            pertrial_depth = struct([]);
            for ii = 1:length(pertrial_mpfr)
                pertrial_depth{ii,1} = pertrial_mpfr{ii,1} - mean(pertrial_bsfr{1,1});
            end

            if isequal(num_targets, 1)
                avg_depth = pertrial_depth{1,1};
            end

            if num_targets >= 2
                if strcmp(tgt_mpfr, 'Max')
                    max_target_center_idx = tgts_per_dir == max(tgts_per_dir);
                    avg_depth = pertrial_depth{max_target_center_idx};
                end
                if strcmp(tgt_mpfr, 'Min')
                    min_target_center_idx = tgts_per_dir == min(tgts_per_dir);
                    avg_depth = pertrial_depth{min_target_center_idx};
                end
                if strcmp(tgt_mpfr, 'Mid')
                    mid_tgt_cntr = tgts_per_dir(ceil(numel(tgts_per_dir)/2));
                    mid_target_center_idx = tgts_per_dir == mid_tgt_cntr;
                    avg_depth = pertrial_depth{mid_target_center_idx};
                end
            end

            %% Put this info in the output arrays
            % Unit names
            all_unit_names{xx,1}(cc,1) = xds.unit_names(uu);

            % Preferred Direction
            units_pref_dir{xx,1}(cc,1) = pref_dir;

            % Target center
            if strcmp(tgt_mpfr, 'Max')
                Tgt_center{xx,1}(cc,1) = max(target_centers(target_dirs == pref_dir));
            end
            if strcmp(tgt_mpfr, 'Min')
                Tgt_center{xx,1}(cc,1) = min(target_centers(target_dirs == pref_dir));
            end
            if strcmp(tgt_mpfr, 'Mid')
                Tgt_center{xx,1}(cc,1) = tgts_per_dir(ceil(numel(tgts_per_dir)/2));
            end

            % Depth of modulation
            depth{xx,1}(cc,1) = mean(avg_depth);

            %% Add to the counter
            cc = cc + 1;
      
        end % End the target direction loop

    end % End the unit loop

    %% Create matrix for all units morning and afternoon

    % Create the excel matrix
    morn_and_noon = cell(length(all_unit_names{xx,1}) + 1, 7);
    % Define the excel matrix headers
    morn_and_noon{1,1} = 'unit_names';
    morn_and_noon{1,2} = 'pref_dir';
    morn_and_noon{1,3} = 'depth';
    morn_and_noon{1,4} = 'tgt_center';

    % Assign the values to the excel matrix
    uu = 2;
    for N = 1:length(all_unit_names{xx,1})
        morn_and_noon{uu,1} = char(all_unit_names{xx,1}(N,1));
        morn_and_noon{uu,2} = units_pref_dir{xx,1}(N,1);
        morn_and_noon{uu,3} = depth{xx,1}(N,1);
        morn_and_noon{uu,4} = Tgt_center{xx,1}(N,1);
        uu = uu + 1;
    end

    %% Save to Excel

    if isequal(Save_Excel, 1)

        % Define the file name
        Monkey_Name = xds.meta.monkey;
        filename = char(strcat(Dates{xx,1}, '_', Monkey_Name, '_', Tasks{xx,1}, '_', Drug_Choice));

        % Save the file
        writecell(morn_and_noon, strcat(Save_Path, filename, '.xlsx'))

    end

end % End the xds loop






