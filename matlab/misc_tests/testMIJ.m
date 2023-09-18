javaaddpath '/Applications/MATLAB_R2023a.app/java/jar/mij.jar'
javaaddpath '/Applications/MATLAB_R2023a.app/java/jar/ij.jar'
MIJ.start;
image = mijread('/Users/lionnt01/Documents/data/feiyue/testSer5ph/batch_1smallImage/Ctrl/8-Ctrl-646,MPM2-488,Ser5ph-Cy3-1001zCorr.tif');

MIJ.closeAllWindows;
MIJ.exit;

% from http://bigwww.epfl.ch/sage/soft/mij/
% MIJ.closeAllWindows	Close all windows
% MIJ.createColor	exports RGB image
% MIJ.createImage	exports other images
% MIJ.error	display a error message in a dialog box
% MIJ.exit	exists MIJ
% MIJ.getCurrentImage	returns a 2D array representing the current image
% MIJ.getImage	returns a 2D array representing the image specified by the title
% MIJ.getColumn	imports column of the ResultsTable
% MIJ.getCurrentTitle	imports the title of the current image
% MIJ.getHistogram	imports the histogram of the current image
% MIJ.getLog	returns the contents of the window log of ImageJ
% MIJ.getListImages	returns the list of opened images
% MIJ.getListColumns	returns the list of columns currently used in the ResultsTable
% MIJ.getResultsTable	imports the ResultsTable
% MIJ.getRoi	imports the current ROI
% MIJ.help	gives a brief description of the MIJ methods
% MIJ.log	print a message in the window console of ImageJ
% MIJ.selectWindow	select a window
% MIJ.setRoi	exports the current ROI
% MIJ.setColumn	exports contents to a column in the ResultsTable
% MIJ.setSlice	set the slice of a stack of image
% MIJ.setThreshold	sets the threshods of the image
% MIJ.setupExt	points to the folder containing ij.jar and plugins and macros folder
% MIJ.showStatus	display a message in the status bar of ImageJ
% MIJ.start	starts MIJ
% MIJ.run	runs command or macro
% MIJ.version	return the MIJ version