The code contained in this folder was written by Rikki Garner to generate the figures in Velle et al. 2023 "A conserved pressure-driven mechanism for regulating cytosolic osmolarity"

To run this code, download the /For_publication/ folder to your computer. Next, add /For_publication/ and all sub-folders to your MATLAB path. To do this, open MATLAB. In the "Current Folder" window panel, navigate to whatever folder you downloaded /For_publication/ into. Right click on /For_publication/ and select Add to Path -> Selected Folders and SubFolders. Next, double-click the /For_publication/ folder to enter that directory. 

This directory contains four main folders: 

1. Figure_Code - The code used to generate the allowable parameters shown in Figure S4A.
	To run this code, double click on /Figure_Code/plot_FigS4A_plotAllowableParameters_CytoPressureModel.m to open this script. Then run the script by clicking "Run" on the topmost panel. The code should plot a range of allowable pore sizes and cytoplasmic pressures, given a set of input parameters specificed in the section "Choose the parameters of the model", and save that plot to Figure_Panels.
2. Figure_Panels - The plots outputted from Figure_Code
3. Data - The raw (input to Data_Analysis_Code) and analyzed (output of Data_Analysis_Code) data files in *.xlsx format.
4. Data_Analysis_Code - The code used to analyze the data (input to Data_Analysis_Code) for generating Fig. 4B-C.
	This code takes in raw traces of vacuole area over time, and determines at each timepoint whether the vacuole is either (a) rapidly expelling its contents or (b) in a paused state with slow or no volume loss. To run this code, open /Data_Analysis_Code/Fig4B_FindRapidLoss_Naegleria.m or /Data_Analysis_Code/Fig4B_FindRapidLoss_Dicty.m and run these scripts by clicking "Run" on the topmost panel. The code should load the raw data, identify paused and unpaused states, plot the results for each vacuole, and then save a new analyzed *.xlsx file. 
5. Accessory_Code - Contains code this author did not write, but are used to run other scripts.
