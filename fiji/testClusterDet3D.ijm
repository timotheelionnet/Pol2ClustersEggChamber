macro ""{
	setBatchMode(true);
	pol2Channel = 4;
	MedianFilterRadiusClusters = 2;
	MedianFilterRadiusBackground = 50;
	maskedNucleoplasmImgTitle = "plasm";
	threshFactor = 3;
	inversion = 1;
	//detectClusters("data","clusterMask",pol2Channel,
	//  		MedianFilterRadiusClusters,MedianFilterRadiusBackground,
	//  		maskedNucleoplasmImgTitle,threshFactor,inversion);
	detectClusters3D("data","clusterMask",pol2Channel,
	  		MedianFilterRadiusClusters,MedianFilterRadiusBackground,
	  		maskedNucleoplasmImgTitle,threshFactor);
	setBatchMode("exit and display");
}

function detectClusters3D(curImg,maskedClustersImgTitle,pol2channel,MedFiltRadiusCluster,MedFiltRadiusBg,maskedNucleoplasmImgTitle,threshF){

	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+pol2Channel);
	run("Median...", "radius="+MedFiltRadiusCluster+" stack");
	rename("pol2median");

	run("Intensity Measurements 2D/3D", "input=pol2median labels="+maskedNucleoplasmImgTitle+" mean stddev");

	// rename the morpholibj output table so ImageJ recognizes it as a Results table
	// the extension of the filename is not transfered to the morpholibJ table name
	Table.rename("pol2median-intensity-measurements", "Results");

	// extract  mean and std
	plasmMean = getResult("Mean",0);
	plasmStd = getResult("StdDev",0);

	// threshold stack
	t = plasmMean + threshF*plasmStd;
	selectWindow("pol2median");
	setThreshold(t, 65535);
	run("Convert to Mask", "method=Default background=Dark black");

}

function detectClusters(curImg,maskedClustersImgTitle,pol2channel,MedFiltRadiusCluster,MedFiltRadiusBg,maskedNucleoplasmImgTitle,threshF){
	verbose = 1;
	// duplicate pol2 channel and median filter to remove salt and pepper noise
	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+pol2Channel);
	run("Median...", "radius="+MedFiltRadiusCluster+" stack");
	rename("pol2median");

	//duplicate nucleoplasm mask
	selectWindow(maskedNucleoplasmImgTitle);
	run("Duplicate...", "duplicate channels=1");
	rename("plasmMask");

	// generate placeholder images 
	getDimensions(w, h, c, nzs, f);
	newImage("clusterMasks", "8-bit Black", w, h, nzs);
	newImage("bgCorrImg", "32-bit Black", w, h, nzs);

	// run through slices 
	for (izs = 0; izs < nzs; izs++) {
	//for (izs = 24; izs < 25; izs++) {
	
		roiManager("reset");
		selectWindow("plasmMask");
		Stack.setPosition(1, izs+1, 1);
		getStatistics(area, mean, min, max, std);

		if (verbose){
			print("slice "+izs+"; minPlasmMask = "+min+"; maxPlasmMask = "+max);
		}
		
		if ((min==0) && (max == 255)){
			// collect avg and std of smoothed pol2 signal
			setThreshold(254, 255);
			run("Create Selection");
			roiManager("Add");
			selectWindow("pol2median");
			Stack.setPosition(1, izs+1, 1);
			roiManager("Select", 0);
			getStatistics(area, mean, min, max, std);
			if (verbose){
				print("slice "+izs+"; area = "+area+"; mean = "+mean+"; std = "+std);
			}
			
			// duplicate smoothed pol2
			run("Select None");
			run("Duplicate...", " ");
			rename("tempPol2");
	
			//duplicate nucleoplasm mask
			selectWindow("plasmMask");
			Stack.setPosition(1, izs+1, 1);
			run("Select None");
			run("Duplicate...", " ");
			run("Divide...", "value=255.000 stack");
			setMinAndMax(0, 65535);
			run("16-bit");
			rename("tempMask");

			// create image that is 0 outside nucleoplasm, smoothed pol2 signal inside
			imageCalculator("Multiply create stack", "tempPol2","tempMask");
			rename("insideNucleoplasm");
	
			// create image that is 0 inside nucleoplasm, <average smoothed pol2 signal within nucleoplasm> outside
			selectWindow("tempMask");
			run("Multiply...", "value="+mean+" stack");
			rename("outsideNucleoplasm");
	
			// sum both images above to generate an image where the holes in the nucleus have the average nucleoplasm intensity
			imageCalculator("Add", "insideNucleoplasm","outsideNucleoplasm");
			rename("filledNucleoplasm");
	
			// median-filter the image just generated and subtract it to remove local background
			selectWindow("filledNucleoplasm");
			run("Duplicate...", "duplicate");
			run("Median...", "radius="+MedFiltRadiusBg);
			rename("medFilteredBackground");
			imageCalculator("Subtract create 32-bit stack", "filledNucleoplasm","medFilteredBackground");
			rename("backgroundSubstracted");
	
			// copy background subtracted image in stack
			run("Select All");
			setPasteMode("Copy");
			run("Copy");
			selectWindow("bgCorrImg");
			Stack.setPosition(1, izs+1, 1);
			setPasteMode("Copy");
			run("Paste");
			
			// compute std of background corrected image
			roiManager("reset")
			selectWindow("plasmMask");
			Stack.setPosition(1, izs+1, 1);
			setThreshold(254, 255);
			run("Create Selection");
			roiManager("Add");
			selectWindow("backgroundSubstracted");
			roiManager("Select", 0);
			
			// set threshold to detect clusters
			getStatistics(area, mean, min, max, std);
			t = threshF*std;
			print("slice "+izs+"; background-corrected std = "+std+ "; threshold = "+t);
			setThreshold(t, 65535);
			run("Convert to Mask", "method=Default background=Dark black");
			setPasteMode("Copy");
			run("Copy");
			selectWindow("clusterMasks");
			Stack.setPosition(1, izs+1, 1);
			setPasteMode("Copy");
			run("Paste");

			// close intermediates
			selectWindow("backgroundSubstracted");
			close();
			selectWindow("tempPol2");
			close();
			selectWindow("filledNucleoplasm");
			close();
			selectWindow("medFilteredBackground");
			close();
			selectWindow("outsideNucleoplasm");
			close();
		}
	}
	selectWindow("clusterMasks");
	// run connected components labeling so that each nucleus gets an individual ID
	run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
	run("glasbey on dark");
	rename(maskedClustersImgTitle);
	
	// close intermediate
	selectWindow("clusterMasks");
	close();
	selectWindow("plasmMask");
	close();
	selectWindow("bgCorrImg");
	close();
	selectWindow("pol2median");
	close();
}