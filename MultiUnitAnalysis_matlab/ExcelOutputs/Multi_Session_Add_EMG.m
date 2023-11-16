function Multi_Session_Add_EMG(Monkey, event, Drug_Choice, Save_Excel)

%% Display the function being used
disp('Multi-Experiment Add EMG Results:');

%% Some of the analysis specifications

% Sorted or unsorted (1 vs 0)
Sorted = 1;

% What muscle groups do you want to look at? ('Flex', 'Ext',
% 'Uln_Dev', 'Rad_Dev', 'Both', 'Grasp', 'Custom', or 'All')
muscle_group = 'Grasp';

% What target direction do you want to use?
dir_choice = 90;

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
    target_center = max(xds_excel.target);

    % Load the xds files
    xds_morn = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Noon');

    % Process the xds files
    Match_The_Targets = 0;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);

    % Zero? (1 = Yes, 0 = No)
    zero_EMG = 1;
    zero_method = 'Prctile';
    EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_group, zero_method, zero_EMG);

    % Normalize EMG? (1 = Yes, 0 = No)
    norm_EMG = 1;
    norm_prctile = 99;
    EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_group, norm_prctile, norm_EMG);

    %% Build the output arrays

    % EMG Indices
    [M] = EMG_Index(xds_morn, muscle_group);

    % EMG names
    EMG_Names = strings;

    % EMG Amplitude
    EMG_amp_morn = NaN(length(unit_names), 1);
    EMG_amp_noon = NaN(length(unit_names), 1);

    % Morning EMG amplitude error
    EMG_amp_err_morn = NaN(length(unit_names), 1);
    % Afternoon EMG amplitude error
    EMG_amp_err_noon = NaN(length(unit_names), 1);

    % EMG amplitude statistics
    EMG_amp_t_test = NaN(length(unit_names), 1);
    EMG_amp_wilcoxon = NaN(length(unit_names), 1);
    EMG_amp_perc = NaN(length(unit_names), 1);
    EMG_amp_cohen_d = NaN(length(unit_names), 1);

    %% Get the peak EMG amplitude
    [~, ~, pertrial_amp_morn] = ...
        EventPeakEMG(xds_morn, muscle_group, event);
    [~, ~, pertrial_amp_noon] = ...
        EventPeakEMG(xds_noon, muscle_group, event);

    %% Extract the target directions & centers
    [target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
    [target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

    phase_idx_morn = intersect(find(target_dirs_morn == dir_choice), find(target_centers_morn == target_center));
    phase_idx_noon = intersect(find(target_dirs_noon == dir_choice), find(target_centers_noon == target_center));

    %% EMG amplitudes at the preferred direction & target center

    pertrial_amp_morn = pertrial_amp_morn(phase_idx_morn);
    pertrial_amp_noon = pertrial_amp_noon(phase_idx_noon);

    %% Loop through the EMG
    for jj = 1:length(M)

        %% Zero the average EMG
        pertrial_amp_morn{1,1}(:,jj) = pertrial_amp_morn{1,1}(:,jj) - EMG_Zero_Factor(jj);
        pertrial_amp_noon{1,1}(:,jj) = pertrial_amp_noon{1,1}(:,jj) - EMG_Zero_Factor(jj);
            
        %% Normalize the average EMG
        pertrial_amp_morn{1,1}(:,jj) = (pertrial_amp_morn{1,1}(:,jj) / EMG_Norm_Factor(jj))*100;
        pertrial_amp_noon{1,1}(:,jj) = (pertrial_amp_noon{1,1}(:,jj) / EMG_Norm_Factor(jj))*100;

        %% Assign the EMG name
        EMG_Names{jj,1} = char(strrep(xds_morn.EMG_names(M(jj)), 'EMG_', ''));

        %% Find the mean & standard error of the EMG amplitude
        EMG_amp_morn(jj,1) = mean(pertrial_amp_morn{1,1}(:,jj));
        EMG_amp_noon(jj,1) = mean(pertrial_amp_noon{1,1}(:,jj));
        EMG_amp_err_morn(jj,1) = std(pertrial_amp_morn{1,1}(:,jj)) / ...
            sqrt(length(pertrial_amp_morn{1,1}(:,jj)));
        EMG_amp_err_noon(jj,1) = std(pertrial_amp_noon{1,1}(:,jj)) / ...
            sqrt(length(pertrial_amp_noon{1,1}(:,jj)));

        %% Check if the EMG amplitude changed significantly

        % EMG amplitude statistics (Unpaired T-Test)
        [~, EMG_amp_t_test(jj,1)] = ttest2(pertrial_amp_morn{1,1}(:,jj), pertrial_amp_noon{1,1}(:,jj));
        % EMG amplitude statistics (Wilcoxon rank sum test)
        [EMG_amp_wilcoxon(jj,1), ~] = ranksum(pertrial_amp_morn{1,1}(:,jj), pertrial_amp_noon{1,1}(:,jj));
        % EMG amplitude percent change
        EMG_amp_perc(jj,1) = (EMG_amp_noon(jj,1) - EMG_amp_morn(jj,1)) / EMG_amp_morn(jj,1);
        % EMG amplitude effect size (Cohen d)
        EMG_amp_cohen_d(jj,1) = Cohen_D(pertrial_amp_morn{1,1}(:,jj), pertrial_amp_noon{1,1}(:,jj));
  
    end % End the EMG amplitude loop

    %% Add the new metric to the excel

    excel_length = length(xds_excel.unit_names);

    % Split the original table
    addition_index = find(strcmp(xds_excel.Properties.VariableNames, 'num_targets'));
    first_half_excel = xds_excel(:,1:addition_index);
    second_half_excel = xds_excel(:,addition_index + 1:end);
    % Create the table additions
    excel_addition = array2table(NaN(excel_length, 9));
    excel_addition.Properties.VariableNames = {'EMG_names', 'EMG_amp_morn', 'EMG_amp_noon', ...
        'EMG_amp_err_morn', 'EMG_amp_err_noon', 'EMG_amp_t_test', 'EMG_amp_wilcoxon', ...
        'EMG_amp_perc', 'EMG_amp_cohen_d'};
    excel_addition.EMG_names = strings(height(excel_addition), 1);
    for ii = 1:length(M)
        excel_addition.EMG_names(ii) = EMG_Names(ii);
    end
    excel_addition.EMG_amp_morn = EMG_amp_morn;
    excel_addition.EMG_amp_noon = EMG_amp_noon;
    excel_addition.EMG_amp_err_morn = EMG_amp_err_morn;
    excel_addition.EMG_amp_err_noon = EMG_amp_err_noon;
    excel_addition.EMG_amp_t_test = EMG_amp_t_test;
    excel_addition.EMG_amp_wilcoxon = EMG_amp_wilcoxon;
    excel_addition.EMG_amp_perc = EMG_amp_perc;
    excel_addition.EMG_amp_cohen_d = EMG_amp_cohen_d;

    % Join the tables
    try
        xds_excel = [first_half_excel excel_addition second_half_excel];
    catch
        xds_excel.EMG_amp_morn = EMG_amp_morn;
        xds_excel.EMG_amp_noon = EMG_amp_noon;
        xds_excel.EMG_amp_err_morn = EMG_amp_err_morn;
        xds_excel.EMG_amp_err_noon = EMG_amp_err_noon;
        xds_excel.EMG_amp_t_test = EMG_amp_t_test;
        xds_excel.EMG_amp_wilcoxon = EMG_amp_wilcoxon;
        xds_excel.EMG_amp_perc = EMG_amp_perc;
        xds_excel.EMG_amp_cohen_d = EMG_amp_cohen_d;
    end


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




