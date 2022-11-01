macro "testSegmentNucleoplasm3D"{
	pol2ChannelNum = 4;
	curImg = getTitle();
	outputImgName = "plasmMask";
	segmentNucleoplasm3D(curImg,pol2ChannelNum,outputImgName);
}

function segmentNucleoplasm3D(curImg,pol2ChannelNum,outputImgName){
	
	selectWindow(curImg);
	getDimensions(w, h, c, nzs, f);
	
	// gnerate a nucleus mask image: threshold nuclei IDs so that all nuclei are white, background is black
	run("Duplicate...", "duplicate channels="+c);
	setThreshold(1, 65535);
	run("Convert to Mask", "method=Default background=Dark black");
	run("16-bit");
	run("Divide...", "value=255.000 stack"); // set white to 1
	rename("masksBinary");

	// duplicate pol2 and multiply it with the [0,1] binary nuclei mask image
	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+pol2ChannelNum);
	rename("pol2tmp");
	imageCalculator("Multiply create stack", "pol2tmp", "masksBinary");
	rename("maskedNucleusIntensity");
	
	// close intermediates
	selectWindow("pol2tmp");
	close();
	selectWindow("masksBinary");
	close();
	
	// make a 0/255 mask for the nucleoplasm alone using autothresholding
	selectWindow("maskedNucleusIntensity");
	//run("Convert to Mask", "method=Default background=Dark black");
	//rename(outputImgName);
	//run("Close-", "stack"); 
}

function segmentNucleoplasm(curImg,pol2ChannelNum,outputImgName){
	
	selectWindow(curImg);
	getDimensions(w, h, c, nzs, f);
	
	// gnerate a nucleus mask image: threshold nuclei IDs so that all nuclei are white, background is black
	run("Duplicate...", "duplicate channels="+c);
	setThreshold(1, 65535);
	run("Convert to Mask", "method=Default background=Dark black");
	run("16-bit");
	run("Divide...", "value=255.000 stack"); // set white to 1
	rename("masksBinary");

	// duplicate pol2 and multiply it with the [0,1] binary nuclei mask image
	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+pol2ChannelNum);
	rename("pol2tmp");
	imageCalculator("Multiply create stack", "pol2tmp", "masksBinary");
	rename("maskedNucleusIntensity");
	
	// close intermediates
	selectWindow("pol2tmp");
	close();
	selectWindow("masksBinary");
	close();
	
	// make a 0/255 mask for the nucleoplasm alone using autothresholding
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
	rename(maskedNucleoplasmImgTitle);
	run("Close-", "stack"); 
}