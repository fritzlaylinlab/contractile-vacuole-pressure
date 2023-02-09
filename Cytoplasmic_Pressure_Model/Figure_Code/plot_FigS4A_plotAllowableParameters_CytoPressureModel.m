%% The following code was written by Rikki Garner to generate the figures 
% in Velle et al. 2023 "A conserved pressure-driven mechanism for
% regulating cytosolic osmolarity"

% The code should plot a range of allowable pore sizes and cytoplasmic pressures, 
% given a set of input parameters specificed in the section 
% "Choose the parameters of the model", and save that plot to Figure_Panels.

%% Clear the systems
close all
clear all

%% Choose the parameters of the model

% Choose the parameters to calculate the resistance to flow through the pore
    % Temperature in Celsius
        T_C = 25;    
    % Boltzman constant in pN nm / K
        k_B = 0.0138;
% Choose the initial shape and size of the vacuole
    % Height of the vacuole in nanometers
    H_V = 10^4;    
    % Initial diameter of the vacuole in nanometers
    D_V_0 = 10^4;
% Choose the time to expel the liquid in ms
    t_E = 10^3;

%% Calculate the derived parameters

% Initial cross-sectional area of the vacuole in nm^2
    A_V_0 = pi*(D_V_0/2)^2;    
% Initial volume of the vacuole in nm^3
    V_V_0 = A_V_0*H_V;
% Calcuate the temperature in K
    T_K = 273.15 + T_C;
% Calculate the dynamic viscosity of water at this temperature in pN ms nm^(-2)
    etaVal = 2.414*10^(-5)*10^(-3)*10^(247.8/(T_K-140));

%% Calculate the corresponding ranges of acceptable parameters

% Choose the pressure across the pore in pN nm^(-2) (= 1 MPa)
    dPVals = [10, 100, 1000]*10^(-6);

% Width of the pore in nanometers
    L_P_Vals = (0:0.01:10)*10^2;

% Preallocate space to store the variables   
    D_P_vals = nan([length(dPVals),length(L_P_Vals)]);

for dPNum = 1:length(dPVals)
    dP = dPVals(dPNum);
    D_P_vals(dPNum,:) = sqrt((4/pi).*sqrt(V_V_0.*8.*pi.*etaVal.*L_P_Vals/(t_E*dP)));
end

% Calculate the areas of the pore in um^2
    A_P_vals = pi*(D_P_vals./2).^2.*10^(-6);   

% Plot the results
    % Choose the figure number
        figure(1)
    % Choose the differnt line stles for each pressure
        plotStyles = {'m-','m--','m:'};
    % Plot the filled region between the min and max pressures
        x2 = [L_P_Vals(2:end), fliplr(L_P_Vals(2:end))];
        inBetween = [D_P_vals(1,(2:end)), fliplr(D_P_vals(end,(2:end)))];
        fill(x2, inBetween, 'black','FaceAlpha',0.2,'HandleVisibility','off');
        hold on;
    % Plot the D_P vs L_P for each pressure choice
        for dPNum = 1:length(dPVals)
        plot(L_P_Vals, D_P_vals(dPNum,:),plotStyles{dPNum},'LineWidth',2)
        end
    % Clean up the plot
        hold off;
        xlabel('Pore depth (nm)')
        ylabel('Pore diameter (nm)')
        legendCell = cellstr(num2str(dPVals'*10^(6), 'P=%2.0f Pa'));
        legend(legendCell,'Location','Northwest')
        title('Allowable parameters given 1 sec dispensing')
       % set(gca,'XScale','log')
       % set(gca,'YScale','log')
        set(gca,'FontName','Helvetica','FontSize',10); 

        % Save the figure 
    % Choose the filename and path for the figure
        destinationTrack = ['C:\Users\Rikki Garner\Documents\Postdoc\'...
            'Research\Collaborations\Velle\Figures\SuppFig6_allowableParameters_cytoPressureModel'];
    % Choose the figure size and position
        figHandle = figure(1);
        figHandle.Position =  [500 400 400 350];
    % Save the file
        exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
    % Save the figure in matlab figure format
        savefig([destinationTrack '.fig'])  
