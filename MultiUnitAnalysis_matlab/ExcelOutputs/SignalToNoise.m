function [SNR_matrix] = SignalToNoise(xds, unit_name, Save_Excel)

%% Find the meta info to load the output excel table

if ~ischar(unit_name)
    % Date
    file_name = xds.meta.rawFileName;
    xtra_info = extractAfter(file_name, '_');
    trial_date = erase(file_name, strcat('_', xtra_info));

    % Task
    if strcmp(xds.meta.task, 'multi_gadget')
        trial_task = 'PG';
    else
        trial_task = 'WS';
    end

    % Monkey
    monkey_name = xds.meta.monkey;

    % Excel file name
    excel_file_name = strcat(trial_date, '_', monkey_name, '_', trial_task);

    % File Path
    file_path = strcat('C:\Users\rhpow\Documents\Grad School\Excel Data\');
    % Files & folders
    dir_path = dir(file_path);
    % Table of only files
    files_in_path = struct2table(dir_path(~([dir_path.isdir])));

    % Selected file
    selec_file = contains(files_in_path.name, excel_file_name);
    output_xds = readtable(strcat(file_path, files_in_path.name{selec_file,1}));

    %% Find the unit of interest

    unit = output_xds.unit_names(unit_name);

    %% Identify the index of the unit
    unit_idx = find(strcmp(xds.unit_names, unit), 1);

elseif ~strcmp(unit_name, 'All')

   unit_idx = find(strcmp(xds.unit_names, unit_name));

end

%% If looping through all units

if strcmp(unit_name, 'All')
    unit_idx = (1:length(xds.spike_waveforms));
end

%% If The Unit Doesn't Exist

if isempty(unit_idx)
    fprintf('%s does not exist \n', unit_name);
    SNR_matrix = NaN;
    return
end

%% Define the exported spreadsheet

SNR_matrix = cell(length(unit_idx)+1 ,4);
SNR_matrix{1,1} = 'Signal:Noise';
SNR_matrix{1,2} = 'Resting Potential RMS';
SNR_matrix{1,3} = 'Average Crest';
SNR_matrix{1,4} = 'Average Trough';

%% Loop through every unit
for N = 1:length(unit_idx)
    % Extracting the waveforms of the designated unit
    waves = xds.spike_waveforms{N};
    
    %% Find the first ten indices of each spike
    rest_potential = zeros(height(waves),10);
    for ii = 1:length(rest_potential)
        rest_potential(ii,:) = waves(ii,1:10);
    end
    
    %% Finding the root mean squared of the noise
    rms_rest_pot = zeros(length(rest_potential),1);
    for ii = 1:length(rest_potential)
        rms_rest_pot(ii) = rms(rest_potential(ii,:));
    end
    
    avg_rms_rest_pot = mean(rms_rest_pot);
    
    %% Finding the mean of the waveforms
    %Calculating the means
    avg_waves = mean(waves);
    
    %% Crests and troughs
    % Amplitude of troughs and crests
    avg_trough = min(avg_waves);
    avg_crest = max(avg_waves);
    % Peak to peak amplitude
    peak_peak_amp = avg_crest - avg_trough;

    SNR = abs(peak_peak_amp / avg_rms_rest_pot);
 
    %% Putting the values in the SNR matrix
    SNR_matrix{N+1, 1} = SNR;
    SNR_matrix{N+1, 2} = avg_rms_rest_pot;
    SNR_matrix{N+1, 3} = avg_crest;
    SNR_matrix{N+1, 4} = avg_trough;
end

%% Save to Excel
if isequal(Save_Excel,  1)
    monkey_name = xds.meta.monkey;
    % Date
    full_trial_date = xds.meta.dateTime;
    xtra_info = extractAfter(full_trial_date, ' ');
    xtra_info = [' ', xtra_info];
    trial_date_with_slash = erase(full_trial_date, xtra_info);
    % Add a zero if the month is single digit
    trial_month = extractBetween(trial_date_with_slash, '/', '/');
    if isequal(length(trial_month), 1)
        double_digit_trial_month = insertBefore(trial_month, trial_month{1,1}, '0');
        trial_date_with_slash = strrep(trial_date_with_slash, trial_month, double_digit_trial_month);
        trial_month = strcat('0', trial_month);
    end
    % Add a zero if the day is single digit
    trial_day = char(extractAfter(trial_date_with_slash, strcat('/', trial_month)));
    if isequal(length(trial_day), 2)
        double_digit_trial_day = insertBefore(trial_day, trial_day(2), '0');
        trial_date_with_slash = strrep(trial_date_with_slash, trial_day, double_digit_trial_day);
    end
    trial_date = erase(trial_date_with_slash, '/');


    filename = char(strcat('Signal2Noise_', monkey_name, '_', trial_date));

    writecell(SNR_matrix, strcat(filename, '.xlsx'))

end





