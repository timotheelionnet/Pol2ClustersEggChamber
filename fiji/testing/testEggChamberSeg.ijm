// macro that generates a segmentation of the 3D space in an eggchamber that does not contain the nuclei for background calculation
// requires as input a channel where each nucleus is segmented and has the pixel value assigned to the ID of the egg chamber it belongs too
macro "testEggChamberSeg"{
	setBatchMode(true);
	testImgName = getTitle();
	outImageName = "segmentedEC";
	segmentEggChamberFromNuclei(testImgName,outImageName);
	setBatchMode("exit and display");
}


function segmentEggChamberFromNuclei(testImgName,outImageName){
	
	// find out how many egg chambers there are
	run("Z Project...", "projection=[Max Intensity]");
	run("Select All");
	getStatistics(area, mean, minImg, nEC, std, histogram);
	print("Found ",nEC," egg chambers to segment.");
	close();
	
	// loop through egg chambers
	for (i = 1; i <= nEC; i++) {
		print("segmenting egg chamber ",i,"...");
		
		// copy stack
		selectWindow(testImgName);
		run("Duplicate...", "duplicate");
		rename("stackCopy"); 
		
		// threshold to keep only the nuclei that belong to the current egg chamber
		setThreshold(i, i);
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=Default background=Light black");
		
		// loop through the 2D slices 
		selectWindow("stackCopy");
		getDimensions(sizeX, sizeY, C, sizeZ, F);
		newImage("conv_"+i, "16-bit", sizeX, sizeY, sizeZ);
		for (j = 0; j < sizeZ; j++) {
			// select current slice of the nuclei stack
			selectWindow("stackCopy");
			Stack.setPosition(1, j+1, 1);
			
			// convexify the nuclei of the current plane
			run("Duplicate...", " ");
			rename("curSlice");
			run("Select All");
			getStatistics(area, mean, minImg, max, std, histogram);
			if(max>0){
				run("Convexify");
			
				// normalize value of the nuclei to match the egg chamber ID
				selectWindow("curSlice-convex");
				run("16-bit");
				run("Divide...", "value=255");
				run("Multiply...", "value="+i);
				
				// add convexified slice to egg chamber stack
				selectWindow("curSlice-convex");
				run("Select All");
				setPasteMode("Copy");
				run("Copy");
	
				selectWindow("conv_"+i);
				Stack.setPosition(1, j+1, 1);
				setPasteMode("Copy");
				run("Paste");	
			}	
			close("curSlice");
			close("curSlice-convex");
		}
		
		// subtract nuclei from convexified stack
		selectWindow("stackCopy");
		run("16-bit");
		run("Divide...", "value=255");
		run("Multiply...", "value="+i);
		imageCalculator("Subtract stack", "conv_"+i,"stackCopy");

		// combine egg chambers into single stack
		if(i>1){
			imageCalculator("Add stack", "conv_1","conv_"+i);
			close("conv_"+i);
		}
		
		// cleanup
		close("stackCopy");

	}
	// set name of output image
	selectWindow("conv_1");
	rename(outImageName);
	
	print("Done segmenting egg chambers.");
	
}
