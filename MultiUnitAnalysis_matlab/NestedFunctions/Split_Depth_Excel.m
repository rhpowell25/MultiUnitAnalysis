function [split_depth_excel, column_names] = Split_Depth_Excel(xds_depth_excel)

%% Build the output array

split_depth_excel = struct([]);
column_names = xds_depth_excel{1}.Properties.VariableNames;

%% Loop through each of experiments
for xx = 1:length(column_names)

    for ii = 1:length(xds_depth_excel)

        split_depth_excel{xx,1}{ii,1} = xds_depth_excel{ii}.(column_names{xx});

    end

end




        






