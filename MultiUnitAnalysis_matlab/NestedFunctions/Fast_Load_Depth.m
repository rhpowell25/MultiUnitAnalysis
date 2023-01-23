function [bsfr_morn, bsfr_noon, depth_morn, depth_noon, bsfr_err_morn, mpfr_err_morn, bsfr_err_noon, mpfr_err_noon, ... 
    depth_t_test, depth_wilcoxon, wave_p_value, nonlin_p_value, fract_contam, pref_dir, num_targets, drug_dose, file_name, all_unit_names] = ...
    Fast_Load_Depth(Drug_Choice, Monkey_Choice, event, Tgt_Choice, Sampling_Params)
%% Import the excel spreadsheets of the selected drug

% Define where the excel spreadsheets are saved
Base_Path = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Excel_Data\', event, '\');
if strcmp(Tgt_Choice, 'Max')
    Data_Path = strcat(Base_Path, 'Max_Targets\');
end
if strcmp(Tgt_Choice, 'Min')
    Data_Path = strcat(Base_Path, 'Min_Targets\');
end
if strcmp(Tgt_Choice, 'All')
    Data_Path = strcat(Base_Path, 'All_Targets\');
end

% Identify all the excel files in the data path
Excel_Path = strcat(Data_Path, '*.xlsx');
Excel_Files = dir(Excel_Path);
% Convert the struc to a table
Excel_Files_In_Path = struct2table(Excel_Files(~([Excel_Files.isdir])));

%% Find the excel files that use the desired drug
if ~strcmp(Drug_Choice, 'All')
    Excel_Drugs = find(contains(Excel_Files_In_Path.name, Drug_Choice));
else
    Excel_Drugs = (1:length(Excel_Files_In_Path.name));
end

%% Find the excel files that are from the desired monkey
if ~strcmp(Monkey_Choice{1,1}, 'All')
    for ii = 1:length(Monkey_Choice)
        if ii == 1
            Excel_Monkey = find(contains(Excel_Files_In_Path.name, Monkey_Choice{ii}));
        else
            Excel_Monkey = cat(1,Excel_Monkey, find(contains(Excel_Files_In_Path.name, Monkey_Choice{ii})));
        end
    end
else
    Excel_Monkey = (1:length(Excel_Files_In_Path.name));
end

%% Find the intersection of Monkey & Drug choices
Excel_Choice = intersect(Excel_Drugs, Excel_Monkey);

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
mpfr_err_morn = struct([]);
% Afternoon error
bsfr_err_noon = struct([]);
mpfr_err_noon = struct([]);

% Depth of modulation statistics (Unpaired T-Test)
depth_t_test = struct([]);

% Depth of modulation statistics (Wilcoxon rank sum test)
depth_wilcoxon = struct([]);

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

% File Name & Dose
drug_dose = strings;
file_name = strings;

%% Loop through each of experiments
for xx = 1:length(Excel_Choice)

    %% Load the output table

    table_path = strcat(Data_Path, Excel_Files_In_Path.name(Excel_Choice(xx)));
    % Continue if the file is open and unsaved
    if contains(Excel_Files_In_Path.name(Excel_Choice(xx)), '~')
        continue
    end

    xds_excel = readtable(char(table_path));

    % Subsample the file according to the paramaters
    [xds_excel] = Subsample_Excel(xds_excel, Sampling_Params);

    %% Put the info in the output arrays

    % Unit names
    all_unit_names{xx,1} = xds_excel.unit_names;

    % Preferred Direction
    pref_dir{xx,1} = xds_excel.pref_dir;

    if ~strcmp(Tgt_Choice, 'All')

        % Morning baseline firing rate
        bsfr_morn{xx,1} = xds_excel.bsfr_morn;
        % Morning baseline firing rate error
        bsfr_err_morn{xx,1} = xds_excel.bsfr_err_morn;
        % Morning depth error
        mpfr_err_morn{xx,1} = xds_excel.mp_err_morn;
        % Morning depth
        depth_morn{xx,1} = xds_excel.depth_morn;
    
        % Afternoon baseline firing rate
        bsfr_noon{xx,1} = xds_excel.bsfr_noon;
        % Afternoon baseline firing rate error
        bsfr_err_noon{xx,1} = xds_excel.bsfr_err_noon;
        % Afternoon depth error
        mpfr_err_noon{xx,1} = xds_excel.mp_err_noon;
        % Afternoon depth
        depth_noon{xx,1} = xds_excel.depth_noon;
    
        % Depth of modulation statistics (Unpaired T-Test)
        depth_t_test{xx,1} = xds_excel.depth_t_test;

        % Depth of modulation statistics (Wilcoxon rank sum test)
        depth_wilcoxon{xx,1} = xds_excel.depth_wilcoxon;
    
        % Wave shape t-test
        wave_p_value{xx,1} = xds_excel.wave_p_value;
        % Nonlinear energy t-test
        nonlin_p_value{xx,1} = xds_excel.nonlin_p_value;

        % ISI Modes
        fract_contam{xx,1} = xds_excel.fract_contam;
    
        % Number of targets
        num_targets{xx,1} = xds_excel.num_targets;
    
        % Drug dose
        drug_dose(xx) = xds_excel.drug_dose_mg_per_kg(1);

    else

        % Morning baseline firing rate
        bsfr_morn{xx,1} = NaN;
        % Morning baseline firing rate error
        bsfr_err_morn{xx,1} = NaN;
        % Morning depth error
        mpfr_err_morn{xx,1} = NaN;
        % Morning depth
        depth_morn{xx,1} = xds_excel.depth;
    
        % Afternoon baseline firing rate
        bsfr_noon{xx,1} = NaN;
        % Afternoon baseline firing rate error
        bsfr_err_noon{xx,1} = NaN;
        % Afternoon depth error
        mpfr_err_noon{xx,1} = NaN;
        % Afternoon depth
        depth_noon{xx,1} = NaN;
    
        % Individual unit t-test
        depth_t_test{xx,1} = NaN;
    
        % Wave shape t-test
        wave_p_value{xx,1} = NaN;
        % Nonlinear energy t-test
        nonlin_p_value{xx,1} = NaN;
        % ISI Modes
        fract_contam{xx,1} = NaN;
    
        % Number of targets
        num_targets{xx,1} = xds_excel.tgt_center;
    
        % Drug dose
        drug_dose(xx) = NaN;
    end

    % File Name
    Excel_File = strrep(char(Excel_Files_In_Path.name(Excel_Choice(xx))), '.xlsx', '');
    file_name(xx,1) = strrep(Excel_File, '_', {' '});
    

end




        






