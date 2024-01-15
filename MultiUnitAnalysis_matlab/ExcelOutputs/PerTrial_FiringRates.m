function PerTrial_FiringRates(Monkey, event, Drug_Choice, Save_Excel)

%% Display the function being used
disp('Per Trial Firing Rates Excel:');

%% Some of the analysis specifications

% Sorted or unsorted (1 vs 0)
Sorted = 1;

Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\');

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug_Choice);

%% Loop through the different experiments
for xx = 1:length(Dates)
    
    % Load the excel file
    [xds_excel] = Load_Excel(Monkey, Dates{xx,1}, Tasks{xx,1});

    nonlin_idxs = xds_excel.nonlin_p_value >= 0.05;
    unit_names = xds_excel.unit_names(nonlin_idxs);

    pref_dirs = xds_excel.pref_dir(nonlin_idxs);
    target_centers = xds_excel.target(nonlin_idxs);

    % Load the xds files
    xds_morn = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Noon');

    % Process the xds files
    Match_The_Targets = 1;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);

    %% Loop through the units  
    for jj = 1:length(unit_names)

        %% Make sure the morning and noon units are identical
        unit_name = char(unit_names(jj));
        pref_dir = pref_dirs(jj);
        target_center = target_centers(jj);

        fprintf('%s \n', unit_name);

        %% Get the firing rates
        % Baseline Firing Rates
        [~, ~, ~, pertrial_bsfr_morn] = BaselineFiringRate(xds_morn, unit_name);
        [~, ~, ~, pertrial_bsfr_noon] = BaselineFiringRate(xds_noon, unit_name);
        % Peak Firing Rates
        [pertrial_mpfr_morn, pertrial_mpfr_noon, ~] = ...
            EventWindow_Morn_v_Noon(xds_morn, xds_noon, unit_name, pref_dir, target_center, event);

        %% Build the output arrays

        if isequal(jj,1)

            % Morning baseline firing rate
            bsfr_morn = zeros(length(pertrial_bsfr_morn{1,1}), length(unit_names));
            % Afternoon baseline firing rate
            bsfr_noon =  zeros(length(pertrial_bsfr_noon{1,1}), length(unit_names));
        
            % Morning peak firing rate
            mpfr_morn = zeros(length(pertrial_mpfr_morn{1,1}), length(unit_names));
            % Afternoon peak firing rate
            mpfr_noon = zeros(length(pertrial_mpfr_noon{1,1}), length(unit_names));

        end

        %% Put the firing rates in the output arrays

        bsfr_morn(:,jj) = pertrial_bsfr_morn{1,1};
        bsfr_noon(:,jj) = pertrial_bsfr_noon{1,1};
        mpfr_morn(:,jj) = pertrial_mpfr_morn{1,1};
        mpfr_noon(:,jj) = pertrial_mpfr_noon{1,1};
  
    end % End the unit loop

    %% Add the new metric to the excel

    excel_length = max(height(bsfr_morn), height(bsfr_noon));
    excel_width = length(unit_names)*4;

    % Create the table
    xds_excel = array2table(NaN(excel_length, excel_width));
    bsfr_morn_vars = string(strcat('bsfr_morn', '_', unit_names'));
    bsfr_noon_vars = string(strcat('bsfr_noon', '_', unit_names'));
    mpfr_morn_vars = string(strcat('mpfr_morn', '_', unit_names'));
    mpfr_noon_vars = string(strcat('mpfr_noon', '_', unit_names'));
    Variable_Names = strcat({bsfr_morn_vars}, {bsfr_noon_vars}, ...
        {mpfr_morn_vars}, {mpfr_noon_vars});
    Variable_Names = string(Variable_Names{1,1});
    xds_excel.Properties.VariableNames = Variable_Names;
    
    % Assign the firing rates to the table
    bsfr_morn_idxs = contains(Variable_Names, 'bsfr_morn');
    xds_excel(1:height(bsfr_morn), bsfr_morn_idxs) = array2table(bsfr_morn);
    bsfr_noon_idxs = contains(Variable_Names, 'bsfr_noon');
    xds_excel(1:height(bsfr_noon), bsfr_noon_idxs) = array2table(bsfr_noon);
    mpfr_morn_idxs = contains(Variable_Names, 'mpfr_morn');
    xds_excel(1:height(mpfr_morn), mpfr_morn_idxs) = array2table(mpfr_morn);
    mpfr_noon_idxs = contains(Variable_Names, 'mpfr_noon');
    xds_excel(1:height(mpfr_noon), mpfr_noon_idxs) = array2table(mpfr_noon);

    %% Save to Excel

    if isequal(Save_Excel, 1)

        % Define the file name
        filename = char(strcat('Per_Trial', Dates{xx,1}, '_', Monkey, '_', ...
            Tasks{xx,1}, '_', Drug_Choice));

        % Save the file
        if ~exist(Save_Path, 'dir')
            mkdir(Save_Path);
        end
        writetable(xds_excel, strcat(Save_Path, filename, '.xlsx'))

    end

end % End the xds loop




