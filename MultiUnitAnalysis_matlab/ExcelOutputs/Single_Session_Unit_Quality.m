function [wave_p_value, nonlin_p_value, fract_contam, all_unit_names] = ...
    Single_Session_Unit_Quality(Monkey, Date, Task, Save_Excel)

%% Display the function being used
disp('Single Experiment Unit Results:');

%% Some of the analysis specifications

% What minimum signal to noise ration? (#, or NaN)
SNR_Minimum = 5;

Save_Path = 'C:\Users\rhpow\Documents\Work\Northwestern\Excel Data\';

%% Build the output arrays

% Wave shape t-test
wave_p_value = struct([]);
% Nonlinear energy t-test
nonlin_p_value = struct([]);
% Fractional Contamination
fract_contam = struct([]);

% Unit Names
all_unit_names = struct([]);

%% Loop through the different experiments
    
% Load the relevant xds file
xds = Load_XDS(Monkey, Date, Task, 1);

% Set the succesful unit counter
cc = 0;

%% Loop through the units  
for jj = 1:length(xds.unit_names)

    %% Skip the unit if the SNR is below the defined minimum

    if ~isnan(SNR_Minimum)
        [SNR_matrix] = SignalToNoise(xds, char(xds.unit_names(jj)), 0);
        SNR = SNR_matrix{2,1};
        if SNR <= SNR_Minimum
            continue
        end
    end

    %% Add to the counter if above the depth minimum & SNR minimum
    cc = cc + 1;

    %% Check if the unit is well sorted
    [~, wave_sort_metric] = WaveShapes(xds, char(xds.unit_names(jj)), 0, 0);
    if isnan(wave_sort_metric)
        cc = cc - 1;
        continue
    end
    [~, nonlin_sort_metric] = NonLinearEnergy(xds, char(xds.unit_names(jj)), 0, 0);
    if isnan(nonlin_sort_metric)
        cc = cc - 1;
        continue
    end
    [Fractional_Contam, ~] = InterstimulusInterval(xds, char(xds.unit_names(jj)), 0, 0);

    %% Put this info in the output arrays
    % Unit names
    all_unit_names{cc,1} = xds.unit_names(jj);

    % Wave shape t-test
    wave_p_value{cc,1} = wave_sort_metric;
    % Nonlinear energy t-test
    nonlin_p_value{cc,1} = nonlin_sort_metric;
    % Fractional Contamination
    fract_contam{cc,1} = Fractional_Contam;

end % End the unit loop

%% Create matrix for all units morning and afternoon

% Create the excel matrix
morn_and_noon = cell(length(all_unit_names) + 1, 4);
% Define the excel matrix headers
morn_and_noon{1,1} = 'unit_names';
morn_and_noon{1,2} = 'wave_p_value';
morn_and_noon{1,3} = 'nonlin_p_value';
morn_and_noon{1,4} = 'fract_contam';

% Assign the values to the excel matrix
uu = 2;
for N = 1:length(all_unit_names)
    morn_and_noon{uu,1} = char(all_unit_names{N,1});
    morn_and_noon{uu,2} = wave_p_value{N,1};
    morn_and_noon{uu,3} = nonlin_p_value{N,1};
    morn_and_noon{uu,4} = fract_contam{N,1};
    uu = uu + 1;
end

%% Save to Excel

if isequal(Save_Excel, 1)

    % Define the file name
    Monkey_Name = xds.meta.monkey;
    filename = char(strcat(Date, '_', Monkey_Name, '_', Task));

    % Save the file
    writecell(morn_and_noon, strcat(Save_Path, filename, '.xlsx'))

end





