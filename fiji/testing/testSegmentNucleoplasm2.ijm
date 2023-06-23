macro "testSegmentNucleoplasm2"{
	curImg = getTitle();
	nucleoplasmSegmentationChannel = 4;
	nucleusMaskChannel = 6;
	outputImgName = "plasmSeg";
	segmentNucleoplasm2(curImg,nucleoplasmSegmentationChannel,nucleusMaskChannel,outputImgName);
	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+nucleoplasmSegmentationChannel);
	rename("curChannel");
	selectWindow("plasmSeg");
	run("16-bit");
	run("Merge Channels...", "c1=curChannel c2=plasmSeg create");
		
}


function segmentNucleoplasm2(inputImg,nucleoplasmSegmentationChannel,nucleusMaskChannel,outputImgName){
	
	run("Duplicate...", "duplicate channels="+nucleoplasmSegmentationChannel);
	run("Median...", "radius=5 stack");
	getDimensions(w, h, c, nzs, f);
	min = 65535;
	max = 0;
	run("Select All");
	for (izs = 0; izs < nzs; izs++) {
		Stack.setPosition(1, izs+1, 1);
		getStatistics(area, mean, curMin, curMax, std);
		if (curMax>max) {
			max = curMax;
		}
		if (curMin<min) {
			min = curMin;
		}
	}
	
	setMinAndMax(min, max);
	run("8-bit");
	run("Auto Local Threshold", "method=Otsu radius=15 parameter_1=0 parameter_2=0 white stack");
	run("Grays");
	run("16-bit");
	rename("plasmMask");
	
	// duplicate the nucleus mask channel
	selectWindow(inputImg);
	run("Duplicate...", "duplicate channels="+nucleusMaskChannel);
	rename("nucMask");
	setThreshold(1, 65535);
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Dark black");
	run("16-bit");
	
	// multiply the nucleus mask by the nucleoplasm mask to erase all signal outside of
	// the nucleus
	imageCalculator("Multiply stack", "plasmMask","nucMask");
	close("nucMask");
	setMinAndMax(0, 255);
	run("8-bit");
	rename(outputImgName);
}