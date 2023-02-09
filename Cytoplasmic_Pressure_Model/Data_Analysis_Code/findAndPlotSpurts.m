%% The following code was written by Rikki Garner to generate the figures 
% in Velle et al. 2023 "A conserved pressure-driven mechanism for
% regulating cytosolic osmolarity"

%% Find the spurts of vacuole area loss

% Calculate the area loss rate
    areaChangeInUM2pS = diff(vacuoleAreaInUM2,1)./diff(timeInSec,1);
% Pull out the time points associated with the area loss rates
    timeForVolumeChangeInS = timeInSec(1:(end-1),:);

% Threshold the area loss rate to find the spurts
    % Choose a threshold in um^2/s
        cutoff = -14;
    % Threshold the area loss rate data to find the spurts
        spurts = (areaChangeInUM2pS<cutoff);

% Pull out only continous spurts that last more than two timesteps 
    % Preallocate space to store continuous spurts
        continuousSpurts = false(size(spurts));
        continuousSpurtIDs = zeros(size(spurts));
    % Loop through each vacuole and find the continuous spurts
    for vacuoleNum = vacuoleID

        % Find continuous strings of unpaused timesteps, and label each with a
        % different number
            [labelImage, n] = bwlabel(spurts(:,vacuoleNum));
        % Pull out the unqiue labels
            labels = unique(labelImage(labelImage~=0));
        % Pull out the number of timesteps for each label
            counts = groupcounts(labelImage(labelImage~=0));
        % Keep only the segments that exist for longer than two timesteps
            tooShortLabels = labels(counts<=2);

        % Loop through each unwanted spurt and delete it
        for labelNum = tooShortLabels'
            % Delete the label
                labelImage(labelImage==labelNum) = 0;    
        end

        % Record the data
            continuousSpurts(find(labelImage>0),vacuoleNum) = 1;
            continuousSpurtIDs(find(labelImage>0),vacuoleNum) = ...
                labelImage(labelImage>0);        

    end

%% Plot the results

    for vacuoleNum = vacuoleID
        % Plot the data for each vacuole separately
            figure(1)
            % Plot the area vs time
                subplot(2,2,1)
                % Plot the raw data
                    plot(timeInSec(:,vacuoleNum),vacuoleAreaInUM2(:,vacuoleNum),'b-o')
                    hold on;
                % Overlay the spurts in a different color
                    for segNum = find(spurts(:,vacuoleNum))'
                        plot(timeInSec(segNum:(segNum+1),vacuoleNum),...
                            vacuoleAreaInUM2(segNum:(segNum+1),vacuoleNum),'r-')
                    end
                % Clean up the plot
                    hold off;
                    ylabel('Area (μm^2)')
                    xlabel('Time (s)')
                    xlim([0 15])
                    ylim([0 100])
            % Plot the area loss rate vs time
                subplot(2,2,2)
                % Plot the raw data
                    plot(timeForVolumeChangeInS(:,vacuoleNum),areaChangeInUM2pS(:,vacuoleNum),'b-o')
                    hold on;
                % Overlay the spurts in a different color
                    plot(timeForVolumeChangeInS(spurts(:,vacuoleNum),vacuoleNum),...
                        areaChangeInUM2pS(spurts(:,vacuoleNum),vacuoleNum),'ro')
                % Overlay the threshold
                    xlim([0 15])
                    xl = get(gca,'Xlim');
                    plot(xl, [cutoff cutoff],'r-')
                % Clean up the plot
                    hold off;
                    ylim([-150 10])
                    xlabel('Time (s)')
                    ylabel('Area loss rate (μm^2/s)')
            % Plot the area vs time, skipping paused timepoints (i.e., for spurts only)
                subplot(2,2,3)
                % Find continuous strings of unpaused timesteps, and label each with a
                % different number
                    [labelImage, n] = bwlabel(continuousSpurts(:,vacuoleNum));
                % Pull out the unqiue labels
                    labels = unique(labelImage(labelImage~=0));
                % Pull out the number of timesteps for each label
                    counts = groupcounts(labelImage(labelImage~=0));
                % Plot the spurts    
                    for labelNum = labels'
                        idx4Area = find(labelImage==labelNum);
                        idx4Area = idx4Area(1):(idx4Area(end)+1);
                        plot(timeInSec((1:(counts(labels==labelNum)+1)),vacuoleNum),...
                            vacuoleAreaInUM2(idx4Area,vacuoleNum),'b-o')
                        hold on;
                    end
                    hold off;
                    ylabel('Area (μm^2)')
                    xlabel('Time (s)')     
            % Plot the area loss rate vs the area (for spurts only)
                subplot(2,2,4)
                % Calculate the mean area between subsequent time points
                    meanArea = (vacuoleAreaInUM2(find(spurts(:,vacuoleNum))+1,...
                        vacuoleNum)+vacuoleAreaInUM2(find(spurts(:,vacuoleNum)),vacuoleNum))/2;
                % Plot the mean area vs the area loss rate
                    plot(meanArea,areaChangeInUM2pS(spurts(:,vacuoleNum)),'b-o')
                % Clean up the plot
                    xlabel('Area (μm^2)')
                    ylabel('Area loss rate (μm^2/s)')
                    drawnow
    
        % Aggregate the data from all particles on a single plot
            figure(2)

            % Plot the area loss rate vs the area (for spurts only)
            subplot(1,2,1)
                % Calculate the mean area between subsequent time points
                    meanArea = (vacuoleAreaInUM2(find(spurts(:,vacuoleNum))+1,...
                        vacuoleNum)+vacuoleAreaInUM2(find(spurts(:,vacuoleNum)),vacuoleNum))/2;
                % Plot the mean area vs the area loss rate
                    plot(meanArea,areaChangeInUM2pS(spurts(:,vacuoleNum)),'-o')
                    hold on;
                % Clean up the plot
                    xlabel('Area (μm^2)')
                    ylabel('Area loss rate (μm^2/s)')
                    drawnow
            % Plot the area vs time, skipping paused timepoints (i.e., for spurts only)
            subplot(1,2,2)
                % Plot the spurts    
                    for labelNum = labels'
                        plot(timeInSec((1:counts(labels==labelNum)),vacuoleNum),...
                            vacuoleAreaInUM2(labelImage==labelNum,vacuoleNum))
                    end
                    ylabel('Area (μm^2)')
                    xlabel('Time (s)')        
                    hold on;
            
                 %   ginput()    

    end

    hold off;

% Plot a histogram of the area reduction rate during spurts

    figure(3)
    histogram(areaChangeInUM2pS(spurts),25)