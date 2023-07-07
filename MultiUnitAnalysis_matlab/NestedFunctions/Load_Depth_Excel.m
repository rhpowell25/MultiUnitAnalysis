function [xds_depth_excel, file_names] = Load_Depth_Excel(Monkey, Sampling_Params)
%% Import the excel spreadsheets of the selected drug

% Define where the excel spreadsheets are saved
Base_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', Sampling_Params.event, '\');
if isequal(Sampling_Params.sorted, 1)
    Data_Path = strcat(Base_Path, 'Sorted\');
elseif isequal(Sampling_Params.sorted, 0)
    Data_Path = strcat(Base_Path, 'Unsorted\');
end


% Identify all the excel files in the data path
Excel_Path = strcat(Data_Path, '*.xlsx');
Excel_Files = dir(Excel_Path);
% Convert the struc to a table
Excel_Files_In_Path = struct2table(Excel_Files(~([Excel_Files.isdir])));

%% Find the excel files that use the desired drug
if ~strcmp(Sampling_Params.drug_choice, 'All')
    Excel_Drugs = find(contains(Excel_Files_In_Path.name, Sampling_Params.drug_choice));
else
    Excel_Drugs = (1:length(Excel_Files_In_Path.name));
end

%% Find the excel files that are from the desired monkey
if ~strcmp(Monkey{1,1}, 'All')
    for ii = 1:length(Monkey)
        if ii == 1
            Excel_Monkey = find(contains(Excel_Files_In_Path.name, Monkey{ii}));
        else
            Excel_Monkey = cat(1,Excel_Monkey, find(contains(Excel_Files_In_Path.name, Monkey{ii})));
        end
    end
else
    Excel_Monkey = (1:length(Excel_Files_In_Path.name));
end

%% Find the intersection of Monkey & Drug choices
Excel_Choice = intersect(Excel_Drugs, Excel_Monkey);

%% Build the output arrays

xds_depth_excel = struct([]);

file_names = strings;

%% Loop through each of experiments

% Initialize the counter
cc = 1;
for xx = 1:length(Excel_Choice)

    %% Load the output table

    if isequal(length(Excel_Choice), 1)
        table_path = strcat(Data_Path, Excel_Files_In_Path.name(Excel_Choice(1)));
        Excel_File = strrep(char(Excel_Files_In_Path.name(Excel_Choice(1))), '.xlsx', '');
    else
        table_path = strcat(Data_Path, Excel_Files_In_Path.name(Excel_Choice(xx)));
        Excel_File = strrep(char(Excel_Files_In_Path.name(Excel_Choice(xx))), '.xlsx', '');
    end
    % Continue if the file is open and unsaved
    if contains(Excel_Files_In_Path.name(Excel_Choice(xx)), '~')
        continue
    end

    temp_depth_excel = readtable(char(table_path));

    % Subsample according to the task
    if isfield(Sampling_Params,'trial_task') && ~strcmp(Sampling_Params.trial_task, 'All')
        if ~contains(Excel_File, Sampling_Params.trial_task)
            continue
        end
    end

    % Subsample according to the number of targets
    if isfield(Sampling_Params,'max_vs_min') && ~isequal(Sampling_Params.max_vs_min, 0)
        if ismember(temp_depth_excel.num_targets, 1)
            continue
        end
    end

    % Subsample the file according to the other paramaters
    [xds_depth_excel{cc}] = Subsample_Excel(temp_depth_excel, Sampling_Params);

    % File Name
    file_names(cc,1) = strrep(Excel_File, '_', {' '});

    % Add to the counter
    cc = cc + 1;

end

%% Merge all the experiments if you selected 'All' for 'trial_sessions')
if isfield(Sampling_Params,'trial_sessions') && strcmp(Sampling_Params.trial_sessions, 'All')

    % Define the merged session table
    merged_session = struct([]);
    merged_variables = xds_depth_excel{1,1}.Properties.VariableNames;
    merged_session{1,1} = array2table(zeros(0, width(merged_variables)));
    merged_session{1,1}.Properties.VariableNames = merged_variables;

    for xx = 1:length(file_names)
        if isempty(xds_depth_excel{xx})
            continue
        end
        if ~iscell(xds_depth_excel{xx}.drug_dose_mg_per_kg(1))
            merged_session{1,1} = [merged_session{1,1}; xds_depth_excel{xx}];
        else
            xds_depth_excel{xx}.drug_dose_mg_per_kg = NaN(height(xds_depth_excel{xx}), 1);
            merged_session{1,1} = [merged_session{1,1}; xds_depth_excel{xx}];
        end
    end

    xds_depth_excel = merged_session;

end


        






