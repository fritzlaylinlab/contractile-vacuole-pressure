#@ File (style="directory") imageFolder
//script to measure all the overlays in a stack
fileList = getFileList(imageFolder); 
for (i = 0; i < lengthOf(fileList); i++) {
    if (endsWith(fileList[i], ".tif")) { 
        open(imageFolder + File.separator + fileList[i]);
        run("Set Scale...", "distance=1 known=1 pixel=1 unit=pixel");
        run("To ROI Manager");
        nRois = roiManager("count");
        for (k = 0; k < nRois; k++) {
        	roiManager("Select", k);
        	roiManager("Measure");
        }
        L = lengthOf(fileList[i]);
        roiSaveName = imageFolder + File.separator + 
                        substring(fileList[i], 0, L - 4) + "_cv.csv";
        saveAs("Results", roiSaveName);
        run("Clear Results");
        close(fileList[i]);
        roiManager("Deselect");
        roiManager("Delete");
    } 
}