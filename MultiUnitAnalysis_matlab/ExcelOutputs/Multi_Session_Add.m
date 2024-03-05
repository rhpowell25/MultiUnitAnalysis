function Multi_Session_Add(Monkey, event, Drug, Save_Excel)

%% Display the function being used
disp('Multi-Experiment Unit Results:');

%% Some of the analysis specifications

Save_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '\Sorted\');

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug);

%% Loop through the different experiments
for xx = 1:length(Dates)

    % Load the excel file
    [xds_excel] = Load_Excel(Monkey, Dates{xx,1}, Tasks{xx,1});

    unit_names = xds_excel.unit_names;
    pref_dirs = xds_excel.pref_dir;
    target_centers = xds_excel.target;

    % Sorted or unsorted (1 vs 0)
    Sorted = 1;

    % Load the xds files
    xds_morn = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx,1}, Tasks{xx,1}, Sorted, 'Noon');

    % Process the xds files
    Match_The_Targets = 0;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);

    %% Build the output arrays

    % Morning alignment times
    alignment_morn = zeros(length(unit_names), 1);
    % Afternoon aligment times
    alignment_noon = zeros(length(unit_names), 1);

    %% Loop through the units  
    for jj = 1:length(unit_names)

        %% Make sure the morning and noon units are identical
        unit_name = char(unit_names(jj));
        pref_dir = pref_dirs(jj);
        target_center = target_centers(jj);

        fprintf('%s \n', unit_name);

        %% Find the alignment times
        [~, alignment_morn(jj)] = EventWindow(xds_morn, unit_name, pref_dir, target_center, event);
        [~, alignment_noon(jj)] = EventWindow(xds_noon, unit_name, pref_dir, target_center, event);

    end

    %% Add the new metric to the excel

    excel_length = length(xds_excel.unit_names);

    % Split the original table
    addition_index = find(strcmp(xds_excel.Properties.VariableNames, 'post_spike_facil'));
    first_half_excel = xds_excel(:,1:addition_index);
    second_half_excel = xds_excel(:,addition_index + 1:end);
    % Create the table additions
    excel_addition = array2table(NaN(excel_length, 2));
    excel_addition.Properties.VariableNames = {'alignment_morn', 'alignment_noon'};
    excel_addition.alignment_morn = alignment_morn;
    excel_addition.alignment_noon = alignment_noon;

    % Join the tables
    xds_excel = [first_half_excel excel_addition second_half_excel];

    %% Save to Excel

    if isequal(Save_Excel, 1)

        % Define the file name
        filename = char(strcat(Dates{xx,1}, '_', Monkey, '_', ...
            Tasks{xx,1}, '_', Drug));

        % Save the file
        if ~exist(Save_Path, 'dir')
            mkdir(Save_Path);
        end
        writetable(xds_excel, strcat(Save_Path, filename, '.xlsx'))

    end

end % End the xds loop




