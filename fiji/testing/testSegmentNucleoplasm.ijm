macro "testSegmentNucleoplasm"{
	run("Close All");
	open("/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/Ctrl/8-Ctrl-646,MPM2-488,Ser5ph-Cy3-1zCorr/nucTIF/nuc62_masks.tif");
	curImg = getTitle();
	nucleoplasmSegmentationChannel = 4;
	outputImgName = "plasmOut";
	segmentNucleoplasm(curImg,nucleoplasmSegmentationChannel,outputImgName);
}



function segmentNucleoplasm(curImg,nucleoplasmSegmentationChannel,outputImgName){
	
	selectWindow(curImg);
	getDimensions(w, h, c, nzs, f);
	
	// threshold nuclei mask channel so that all nuclei are white, background is black
	run("Duplicate...", "duplicate channels="+c);
	setThreshold(1, 65535);
	run("Convert to Mask", "method=Default background=Dark black");
	run("16-bit");
	run("Divide...", "value=255.000 stack"); // set white to 1
	rename("masksBinary");

	// duplicate signal channel and multiply it with the [0,1] binary nuclei mask image
	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+nucleoplasmSegmentationChannel);
	rename("pol2tmp");
	imageCalculator("Multiply create stack", "pol2tmp", "masksBinary");
	rename("maskedNucleusIntensity");
	
	// close intermediates
	close("pol2tmp");
	close("masksBinary");
	
	// make a 0/255 mask for the nucleoplasm alone using autothresholding
	// of the signal channel only within the nucleus mask
	selectWindow("maskedNucleusIntensity");
	getDimensions(w, h, c, nzs, f);
	
	for (izs = 0; izs < nzs; izs++) {
		Stack.setPosition(1, izs+1, 1);
		getStatistics(area, mean, min, max, std);
		if (max!=0) {
			run("Auto Threshold", "method=Default ignore_black white");
		}
	}
	run("8-bit");
	rename(outputImgName);
	//%run("Close-", "stack"); 
}