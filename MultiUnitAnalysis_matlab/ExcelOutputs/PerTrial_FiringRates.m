function PerTrial_FiringRates(Monkey, event, Drug, Save_Excel)

%% Display the function being used
disp('Per Trial Firing Rates Excel:');

%% Some of the analysis specifications

% Sorted or unsorted (1 vs 0)
Sorted = 1;

Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\');

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug);

%% Loop through the different experiments
for xx = 1:length(Dates)
    
    % Load the excel file
    [xds_excel] = Load_Excel(Monkey, Dates{xx,1}, Tasks{xx,1});

    nonlin_idxs = xds_excel.nonlin_p_val >= 0.05;
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

            % Morning depth of modulation
            depth_morn = zeros(length(pertrial_mpfr_morn{1,1}), length(unit_names));
            % Afternoon depth of modulation
            depth_noon = zeros(length(pertrial_mpfr_noon{1,1}), length(unit_names));

        end

        %% Put the firing rates in the output arrays

        depth_morn(:,jj) = pertrial_mpfr_morn{1,1} - mean(pertrial_bsfr_morn{1,1});
        depth_noon(:,jj) = pertrial_mpfr_noon{1,1} - mean(pertrial_bsfr_noon{1,1});
  
    end % End the unit loop

    %% Add the new metric to the excel

    excel_length = max(height(depth_morn), height(depth_noon));
    excel_width = length(unit_names)*2;

    % Create the table
    xds_excel = array2table(NaN(excel_length, excel_width));
    depth_morn_vars = string(strcat('depth_morn', '_', unit_names'));
    depth_noon_vars = string(strcat('depth_noon', '_', unit_names'));
    Variable_Names = strcat({depth_morn_vars}, {depth_noon_vars});
    Variable_Names = string(Variable_Names{1,1});
    xds_excel.Properties.VariableNames = Variable_Names;
    
    % Assign the firing rates to the table
    depth_morn_idxs = contains(Variable_Names, 'depth_morn');
    xds_excel(1:height(depth_morn), depth_morn_idxs) = array2table(depth_morn);
    depth_noon_idxs = contains(Variable_Names, 'depth_noon');
    xds_excel(1:height(depth_noon), depth_noon_idxs) = array2table(depth_noon);

    %% Save to Excel

    if isequal(Save_Excel, 1)

        % Define the file name
        filename = char(strcat('PerTrial_', Dates{xx,1}, '_', Monkey, '_', ...
            Tasks{xx,1}, '_', Drug));

        % Save the file
        if ~exist(Save_Path, 'dir')
            mkdir(Save_Path);
        end
        writetable(xds_excel, strcat(Save_Path, filename, '.xlsx'))

    end

end % End the xds loop




