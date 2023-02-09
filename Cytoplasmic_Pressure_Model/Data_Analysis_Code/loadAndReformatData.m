%% The following code was written by Rikki Garner to generate the figures 
% in Velle et al. 2023 "A conserved pressure-driven mechanism for
% regulating cytosolic osmolarity"
    
% Load the data
    vacuoleEmptyingData = readtable(dataFilePath,'Sheet', sheetName);

% Reformat the data
    % Create a numTPxnumVacuoles array containing the time points in seconds
        timeInSec = table2array(vacuoleEmptyingData(:,2:3:end));
    % Create a numTPxnumVacuoles array containing the vacuole area for
    % each vacuole and time point
        vacuoleAreaInUM2 = table2array(vacuoleEmptyingData(:,1:3:end));
    % Create identifiers for each vacuole
        pumpID = 1:size(vacuoleAreaInUM2,2);
        vacuoleID = 1:size(vacuoleAreaInUM2,2);
        experimentID = ones(size(pumpID));