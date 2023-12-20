function Multi_Session_Add_Phase(Monkey, event, Drug_Choice, Save_Excel)

%% Display the function being used
disp('Multi-Experiment Add Phase Results:');

%% Some of the analysis specifications

% Sorted or unsorted (1 vs 0)
Sorted = 1;

if isequal(Sorted, 1)
    Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '\Sorted\');
else
    Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '\Unsorted\');
end

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug_Choice);

%% Loop through the different experiments
for xx = 1:length(Dates)
    
    % Load the excel file
    [xds_excel] = Load_Excel(Monkey, Dates{xx,1}, Tasks{xx,1});

    unit_names = xds_excel.unit_names;
    pref_dirs = xds_excel.pref_dir;
    target_centers = xds_excel.target;

    % Load the xds files
    xds_morn = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Noon');

    % Process the xds files
    Match_The_Targets = 1;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);

    %% Build the output arrays

    % Morning ramp firing rate
    ramp_morn = zeros(length(unit_names), 1);
    % Afternoon ramp firing rate
    ramp_noon = zeros(length(unit_names), 1);
    % Morning TgtHold firing rate
    TgtHold_morn = zeros(length(unit_names), 1);
    % Afternoon TgtHold firing rate
    TgtHold_noon = zeros(length(unit_names), 1);

    % Morning ramp firing rate error
    ramp_err_morn = zeros(length(unit_names), 1);
    % Morning TgtHold firing rate error
    TgtHold_err_morn = zeros(length(unit_names), 1);
    % Afternoon ramp firing rate error
    ramp_err_noon = zeros(length(unit_names), 1);
    % Afternoon TgtHold firing rate error
    TgtHold_err_noon = zeros(length(unit_names), 1);

    % Ramp phase statistics
    ramp_p_value_t_test = zeros(length(unit_names), 1);
    ramp_p_value_wilcoxon = zeros(length(unit_names), 1);
    ramp_effect_perc = zeros(length(unit_names), 1);
    ramp_effect_cohen_d = zeros(length(unit_names), 1);

    % TgtHold phase statistics
    TgtHold_p_value_t_test = zeros(length(unit_names), 1);
    TgtHold_p_value_wilcoxon = zeros(length(unit_names), 1);
    TgtHold_effect_perc = zeros(length(unit_names), 1);
    TgtHold_effect_cohen_d = zeros(length(unit_names), 1);

    %% Loop through the units  
    for jj = 1:length(unit_names)

        %% Make sure the morning and noon units are identical
        unit_name = char(unit_names(jj));
        pref_dir = pref_dirs(jj);
        target_center = target_centers(jj);

        fprintf('%s \n', unit_name);

        %% Get the ramp phase firing rates
        [avg_ramp_morn, ~, err_ramp_morn, pertrial_ramp_morn] = ... 
            RampFiringRate(xds_morn, unit_name);
        [avg_ramp_noon, ~, err_ramp_noon, pertrial_ramp_noon] = ... 
            RampFiringRate(xds_noon, unit_name);

        %% Get the target hold firing rates
        [avg_TgtHold_morn, ~, ~, pertrial_TgtHold_morn] = ...
            EventPeakFiringRate(xds_morn, unit_name, 'TgtHold');
        [avg_TgtHold_noon, ~, ~, pertrial_TgtHold_noon] = ...
            EventPeakFiringRate(xds_noon, unit_name, 'TgtHold');

        %% Find the standard error of the target hold firing rates

        err_TgtHold_morn = zeros(length(pertrial_TgtHold_morn),1);
        err_TgtHold_noon = zeros(length(pertrial_TgtHold_noon),1);
        for ii = 1:length(avg_TgtHold_morn)
            err_TgtHold_morn(ii,1) = std(pertrial_TgtHold_morn{ii,1}) / ...
                sqrt(length(pertrial_TgtHold_morn{ii,1}));
            err_TgtHold_noon(ii,1) = std(pertrial_TgtHold_noon{ii,1}) / ...
                sqrt(length(pertrial_TgtHold_noon{ii,1}));
        end

        %% Extract the target directions & centers
        [target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
        [target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

        phase_idx_morn = intersect(find(target_dirs_morn == pref_dir), ...
            find(target_centers_morn == target_center));
        phase_idx_noon = intersect(find(target_dirs_noon == pref_dir), ...
            find(target_centers_noon == target_center));

        %% Phasic firing rates at the preferred direction & target center

        avg_ramp_morn = avg_ramp_morn(phase_idx_morn);
        avg_ramp_noon = avg_ramp_noon(phase_idx_noon);

        avg_TgtHold_morn = avg_TgtHold_morn(phase_idx_morn);
        avg_TgtHold_noon = avg_TgtHold_noon(phase_idx_noon);

        pertrial_ramp_morn = pertrial_ramp_morn(phase_idx_morn);
        pertrial_ramp_noon = pertrial_ramp_noon(phase_idx_noon);

        pertrial_TgtHold_morn = pertrial_TgtHold_morn(phase_idx_morn);
        pertrial_TgtHold_noon = pertrial_TgtHold_noon(phase_idx_noon);
        
        %% Standard error of the phasic firing rates at the preferred direction & target center

        err_ramp_morn = err_ramp_morn(phase_idx_morn);
        err_ramp_noon = err_ramp_noon(phase_idx_noon);

        err_TgtHold_morn = err_TgtHold_morn(phase_idx_morn);
        err_TgtHold_noon = err_TgtHold_noon(phase_idx_noon);

        %% Check if the unit's depth of modulation changed significantly

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

        %% Put this info in the output arrays

        % Morning ramp firing rate
        ramp_morn(jj,1) = avg_ramp_morn;
        % Afternoon ramp firing rate
        ramp_noon(jj,1) = avg_ramp_noon;
        % Morning TgtHold firing rate
        TgtHold_morn(jj,1) = avg_TgtHold_morn;
        % Afternoon TgtHold firing rate
        TgtHold_noon(jj,1) = avg_TgtHold_noon;

        % Morning ramp firing rate error
        ramp_err_morn(jj,1) = err_ramp_morn;
        % Morning TgtHold firing rate error
        TgtHold_err_morn(jj,1) = err_TgtHold_morn;
        % Afternoon ramp firing rate error
        ramp_err_noon(jj,1) = err_ramp_noon;
        % Afternoon TgtHold firing rate error
        TgtHold_err_noon(jj,1) = err_TgtHold_noon;

        % Ramp phase statistics
        ramp_p_value_t_test(jj,1) = ramp_t_test;
        ramp_p_value_wilcoxon(jj,1) = ramp_wilcoxon;
        ramp_effect_perc(jj,1) = ramp_perc;
        ramp_effect_cohen_d(jj,1) = ramp_cohen_d;

        % TgtHold phase statistics
        TgtHold_p_value_t_test(jj,1) = TgtHold_t_test;
        TgtHold_p_value_wilcoxon(jj,1) = TgtHold_wilcoxon;
        TgtHold_effect_perc(jj,1) = TgtHold_perc;
        TgtHold_effect_cohen_d(jj,1) = TgtHold_cohen_d;
  
    end % End the unit loop

    %% Add the new metric to the excel

    excel_length = length(xds_excel.unit_names);

    % Split the original table
    addition_index = find(strcmp(xds_excel.Properties.VariableNames, 'pref_dir'));
    first_half_excel = xds_excel(:,1:addition_index);
    second_half_excel = xds_excel(:,addition_index + 1:end);
    % Create the table additions
    excel_addition = array2table(NaN(excel_length, 16));
    excel_addition.Properties.VariableNames = {'ramp_morn', 'ramp_noon', ...
        'ramp_err_morn', 'ramp_err_noon', 'TgtHold_morn', 'TgtHold_noon', ...
        'TgtHold_err_morn', 'TgtHold_err_noon', 'ramp_t_test', 'ramp_wilcoxon', ...
        'ramp_perc', 'ramp_cohen_d', 'TgtHold_t_test', 'TgtHold_wilcoxon', ...
        'TgtHold_perc', 'TgtHold_cohen_d'};
    excel_addition.ramp_morn = ramp_morn;
    excel_addition.ramp_noon = ramp_noon;
    excel_addition.ramp_err_morn = ramp_err_morn;
    excel_addition.ramp_err_noon = ramp_err_noon;
    excel_addition.TgtHold_morn = TgtHold_morn;
    excel_addition.TgtHold_noon = TgtHold_noon;
    excel_addition.TgtHold_err_morn = TgtHold_err_morn;
    excel_addition.TgtHold_err_noon = TgtHold_err_noon;
    excel_addition.ramp_t_test = ramp_p_value_t_test;
    excel_addition.ramp_wilcoxon = ramp_p_value_wilcoxon;
    excel_addition.ramp_perc = ramp_effect_perc;
    excel_addition.ramp_cohen_d = ramp_effect_cohen_d;
    excel_addition.TgtHold_t_test = TgtHold_p_value_t_test;
    excel_addition.TgtHold_wilcoxon = TgtHold_p_value_wilcoxon;
    excel_addition.TgtHold_perc = TgtHold_effect_perc;
    excel_addition.TgtHold_cohen_d = TgtHold_effect_cohen_d;

    % Join the tables
    xds_excel = [first_half_excel excel_addition second_half_excel];

    %% Save to Excel

    if isequal(Save_Excel, 1)

        % Define the file name
        filename = char(strcat(Dates{xx,1}, '_', Monkey, '_', ...
            Tasks{xx,1}, '_', Drug_Choice));

        % Save the file
        if ~exist(Save_Path, 'dir')
            mkdir(Save_Path);
        end
        writetable(xds_excel, strcat(Save_Path, filename, '.xlsx'))

    end

end % End the xds loop




