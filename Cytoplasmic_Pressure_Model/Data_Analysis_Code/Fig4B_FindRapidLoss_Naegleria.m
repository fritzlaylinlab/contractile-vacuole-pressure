%% The following code was written by Rikki Garner to generate the figures 
% in Velle et al. 2023 "A conserved pressure-driven mechanism for
% regulating cytosolic osmolarity"

% The code should load the raw data, identify paused and unpaused states, 
% plot the results for each vacuole, and then save a new analyzed *.xlsx file. 

% Clear the system
    clear all;
    close all;

% Choose the file path to the data
    dataFilePath = ['\Data\Values for modeling.xlsx'];
% Choose the sheet name in the file
    sheetName = 'emptying faster frame rate';

% Load and reformat the data
    loadAndReformatData    

% Find the moments of rapid vacuole area loss
    findAndPlotSpurts

% Export the data
    % Choose the sheet name
        newSheetName = 'EFFR - ResetTime';
    % Save the results
        saveResults