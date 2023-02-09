%% The following code was written by Rikki Garner to generate the figures 
% in Velle et al. 2023 "A conserved pressure-driven mechanism for
% regulating cytosolic osmolarity"

%% Choose the file path to the data
    % Choose a file name and path
        % Find the file prefix
            % Find the last index before the filetype extension
                lastIdx = find(dataFilePath=='.',1,'last')-1;
            % Pull out the file prexi
                filePathPrefix = dataFilePath(1:lastIdx);
        % Add an addendum and filetype extension
            saveFilePath = [filePathPrefix '_Analyzed.xlsx'];
    % If the file doesn't exist yet, copy the original file to a file with
    % this name
        if ~isfile(saveFilePath)
            status = copyfile(dataFilePath,saveFilePath);
        end     

% Add spurt boolean (1 is spurt, 0 is paused) to the original spreadsheet
    
    % Pull out the associated column numbers for each data point
        colNums = vacuoleID*3;
    % Loop through each vacuole and store the data
    for vacuoleNum = vacuoleID
        nonNan = ~isnan(areaChangeInUM2pS(:,vacuoleNum));
        unPaused = spurts(:,vacuoleNum);
        xlswrite(saveFilePath,unPaused(nonNan),sheetName,[xlscol(colNums(vacuoleNum)) '3'])
        xlswrite(saveFilePath,{'Unpaused?'},sheetName,[xlscol(colNums(vacuoleNum)) '2'])
    end


%% Create a new spreadsheet with the area vs time, skipping paused timepoints (i.e., for spurts only)

    % Initialize the column numbers
        colNum = 1;

    for vacuoleNum = vacuoleID
    

        % Find continuous strings of unpaused timesteps, and label each with a
        % different number
            [labelImage, n] = bwlabel(continuousSpurts(:,vacuoleNum));
        % Pull out the unqiue labels
            labels = unique(labelImage(labelImage~=0));
        % Pull out the number of timesteps for each label
            counts = groupcounts(labelImage(labelImage~=0));
    
        % Label the first row with the vacuole number
            xlswrite(saveFilePath,{sprintf('Vacuole %1.0f',vacuoleNum)},...
                newSheetName,[xlscol(colNum) '1'])
    
        % Initialize the number of spurts
            spurtNum = 0;
    
        % For each label
        for labelNum = labels'
    
            % Update the number of spurts
                spurtNum = spurtNum+1;
    
            % Add unpaused region number label 
                xlswrite(saveFilePath,{sprintf('Spurt %1.0f',spurtNum)},newSheetName,[xlscol(colNum) '2'])
    
            % Add the area vector
                % Pull out the indices to plot
                    idx4Area = find(labelImage==labelNum);
                    idx4Area = idx4Area(1):(idx4Area(end)+1);
                % Add the column label
                    xlswrite(saveFilePath,{'Area (um^2)'},newSheetName,[xlscol(colNum) '3'])
                % Add the data
                    xlswrite(saveFilePath,vacuoleAreaInUM2(idx4Area,vacuoleNum),newSheetName,[xlscol(colNum) '4'])

            % Add the time vector
                % Update the column counter
                    colNum = colNum+1;
                % Add the time column     
                    xlswrite(saveFilePath,{'Reset Time (sec)'},newSheetName,[xlscol(colNum) '3'])
                    xlswrite(saveFilePath,timeInSec((1:(counts(labels==labelNum)+1)),vacuoleNum),newSheetName,[xlscol(colNum) '4'])
    
            
            % Update the column counter
                colNum = colNum+2;
    
        end


    end