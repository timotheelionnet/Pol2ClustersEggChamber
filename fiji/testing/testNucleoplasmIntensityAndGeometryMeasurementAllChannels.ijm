macro "testNucleoplasmIntensityAndGeometryMeasurementAllChannels"{
	run("Close All");
	imgPath = "/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/TRI/9-TRI-646,MPM2-488,Ser5ph-Cy3-1zCorr/nucTIF/nuc4_masks.tif";
	open(imgPath);
	imgOutName = getTitle();
	csvSaveDir = "/Users/lionnt01/Dropbox/data/junk/";
	maskChannel = 6;
	label = "4";
	nucleoplasmIntensityAndGeometryMeasurementAllChannels(
		csvSaveDir,"nuc"+label,imgOutName,maskChannel,maskChannel+1);
}

// takes an image (inputWindowName), duplicates the channel with nucleoplasm mask, 
// then measures 1) intensity etc and 2) geometry on all Channels
// and saves the resulting tables in the folder savePath using names like 
// C1_rootName_plasmInt.csv, C2_rootName_plasmInt.csv, etc for intensity features
// single file rootName-nucleoplasmGeom.csv for geometry features 
// generates a mask for nucleoli (inverse of nucleoplasm within nucleus mask)
// and perform the same measurements of geometry and intensity for the nucleoli mask.
function nucleoplasmIntensityAndGeometryMeasurementAllChannels(savePath,rootName,inputWindowName,
												nucChannel,plasmChannel){
	selectWindow(inputWindowName);
	
	// duplicate the nucleoplasm mask channel and perform geometry measurements on it.
	run("Duplicate...", "duplicate channels="+plasmChannel);
	plasmDuplicate = "tmpPlasmDuplicate";
	rename(plasmDuplicate);
	run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity"+
			" euler_number bounding_box centroid equivalent_ellipsoid ellipsoid_elongations"+
			" max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	
	// save plasm geometry in csv file
	saveAs("Results",savePath+rootName+"_plasmGeom.csv");
	run("Close");
	
	// generate nucleoli mask - i.e. inverse of nucleoplasm within nucleus
	selectWindow(inputWindowName);
	run("Duplicate...", "duplicate channels="+nucChannel);
	nucDuplicate = "tmpNucDuplicate";
	rename(nucDuplicate);
	normalizePixelValuesToBitDepth(nucDuplicate);
	selectWindow(plasmDuplicate);
	run("Duplicate...", "duplicate");
	t = getTitle();
	normalizePixelValuesToBitDepth(t);
	run("Invert", "stack"); // invert nucleoplasm mask
	run("Divide...", "value=65535 stack"); // need to divide because otherwise the image mulitplication later on
	// gives rise to 65,535 * 65,535 which Fiji interprets as zero.
	rename("plasmInv");
	imageCalculator("Multiply stack", "plasmInv",nucDuplicate); // multiply inverse nucleoplasm and nucleus masks
	nucleoliMask = "nucleoliMask";
	rename(nucleoliMask);
	if(isOpen(nucDuplicate)){
		close(nucDuplicate);
	}
	if(isOpen("plasmInv")){
		close("plasmInv");
	}
	//set nucleoli mask value to the same value as nucleoplasm mask
	selectWindow(nucleoliMask);
	run("Z Project...", "projection=[Max Intensity]");
	rename("tmpMax");		
	getStatistics(area, mean, min, maxNucleoli, std, histogram);
	close("tmpMax");
	
	selectWindow(plasmDuplicate);
	run("Z Project...", "projection=[Max Intensity]");
	rename("tmpMax");		
	getStatistics(area, mean, min, maxPlasm, std, histogram);
	close("tmpMax");
	
	selectWindow(nucleoliMask);
	run("Divide...", "value="+maxNucleoli+" stack");	
	run("Multiply...", "value="+maxPlasm+" stack");	
	
	// measure geometry of nucleoli
	run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity"+
			" euler_number bounding_box centroid equivalent_ellipsoid ellipsoid_elongations"+
			" max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	
	// save nucleoli geometry in csv file
	saveAs("Results",savePath+rootName+"_nucleoliGeom.csv");
	run("Close");
	
	// loop through channels and collect intensity measurements over nucleoplasm and nucleoli
	selectWindow(inputWindowName);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	for (i = 1; i <= C; i++) {
	    selectWindow(inputWindowName);
	    run("Duplicate...", "duplicate channels="+i);
	    curChannelDuplicate = "tmpChannelDuplicate";
		rename(curChannelDuplicate);
		
		// nucleoplasm measurements
		run("Intensity Measurements 2D/3D", "input="+curChannelDuplicate+" labels="+plasmDuplicate+ 
			" mean stddev max min median mode skewness kurtosis numberofvoxels volume"+
			" neighborsmean neighborsstddev neighborsmax neighborsmin neighborsmedian"+
			" neighborsmode neighborsskewness neighborskurtosis");
		
		saveAs(curChannelDuplicate+"-intensity-measurements",
					savePath+"C"+i+"_"+rootName+"_plasmInt.csv");
		run("Close");
		
		// nucleoli measurements
		run("Intensity Measurements 2D/3D", "input="+curChannelDuplicate+" labels="+nucleoliMask+ 
			" mean stddev max min median mode skewness kurtosis numberofvoxels volume"+
			" neighborsmean neighborsstddev neighborsmax neighborsmin neighborsmedian"+
			" neighborsmode neighborsskewness neighborskurtosis");
		
		saveAs(curChannelDuplicate+"-intensity-measurements",
					savePath+"C"+i+"_"+rootName+"_nucleoliInt.csv");
		run("Close");
		
		close(curChannelDuplicate);
	}
	
	print("Nucleoplasm Intensity and Geometry results saved");
	close(plasmDuplicate);
	close(nucleoliMask);
}

// takes each channel of the image and applies an affine transformation 
// so that the min value of the stack (or 2D image) in each channel is zero and the max
// value is the bit depth (255 for 8-bit and 65535 for 16 bit).
// replaces the original image with the renormalized one.
function normalizePixelValuesToBitDepth(imgName){
	selectWindow(imgName);
	bd = bitDepth();
	if((bd!= 8) && (bd!=16)){
		print("img "+imgName+" is neither 8 nor 16 bit, cannot normalize.");
		return;
	}
	// set range of renormalized image to entire bitdepth
	bitMax = pow(2, bd)-1;
	
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	if(C>1){
		print("image has "+C+" channels. will normalize each channel to bit depth.");
	}
	for (i = 1; i <= C; i++) {
		selectWindow(imgName);
		run("Duplicate...", "duplicate channels="+i);
		rename("tmpChannel");
		
		// find the min and max of the image
		if(sizeZ>1){
			run("Z Project...", "projection=[Max Intensity]");
			getStatistics(area, mean, min, imgMax, std);
			close();
			run("Z Project...", "projection=[Min Intensity]");
			getStatistics(area, mean, imgMin, max, std);
			close();
		}else {
			getStatistics(area, mean, imgMin, imgMax, std);
		}
		selectWindow("tmpChannel");
		// set min pixel(s) to zero
		run("Subtract...", "value="+imgMin+" stack");
		
		// normalize to bit depth
		if(imgMin < imgMax){
			multFactor = bitMax/(imgMax-imgMin);
			if(sizeZ>1){
				run("Multiply...", "value="+multFactor+" stack");	
			}else {
				run("Multiply...", "value="+multFactor);
			}
		}
		
		if(i==1){
			selectWindow("tmpChannel");
			rename("renormImg");
		}else{
			addChannelToImg("renormImg","tmpChannel","renormImg",0);
		}
	}
	
	// replace original image with new one
	close(imgName);
	selectWindow("renormImg");
	rename(imgName);
}
