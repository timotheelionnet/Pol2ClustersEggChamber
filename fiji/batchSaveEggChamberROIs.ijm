// this macro runs through a folder of raw data and loads each stack sequentially,
// prompting the user to segment egg chambers manually. This is done using the free
// drawing tool to generate ROIs and then enter them for each stack into the ROI manager.
// the macro then saves the ROIs as 2D masks in dedicated subfolders for each image.

// tracing tips: use the plane where the egg chamber is the widest and trace the boundaries. 
// it's ok if egg chambers overlap a bit. The idea is that a nuclei will be assigned to an
// egg chamber if > 90 % of its volume falls within the 2D mask (extended vertically).

macro "batchSaveEggChamberROIs"{
	run("Close All");
	// where data comes from, where segmentations go to.
	inFolder = getDirectory("choose the input directory");
	outFolder = getDirectory("choose the output directory");
	
	// collect names of image files in the input dir and its subfolders and
	// store into outSubDirList & fileList, so that the path to file i is
	// inFolder+outSubDirList[i]+fileList[i].
	dirList = getFileList(inFolder);
	fileList = newArray(1000);
	outSubDirList = newArray(1000);
	ctr = 0;
	print(" ");
	for (i = 0; i < lengthOf(dirList); i++) {
	    if (endsWith(dirList[i], ".tif")) { 
	        fileList[ctr] = dirList[i];
	        outSubDirList[ctr] = "";
	        print("fname: "+fileList[ctr]+"; dir: "+outSubDirList[ctr]+"; ctr: "+ctr);
	        ctr = ctr+1;
	    } else{
	    	subDirList = getFileList(inFolder+dirList[i]);
	    	if (subDirList.length>0){
	    		for (j = 0; j < lengthOf(subDirList); j++) {
	    			if (endsWith(subDirList[j], ".tif")) { 
				        fileList[ctr] = subDirList[j];
				        outSubDirList[ctr] = dirList[i];
				        print("fname: "+fileList[ctr]+"; dir: "+outSubDirList[ctr]+"; ctr: "+ctr);
				        ctr = ctr+1;
				    } 
	    		}
	    	}
	    }
	}
	fileList = Array.trim(fileList, ctr);
	outSubDirList = Array.trim(outSubDirList, ctr);
	
	// create output subfolders if needed, i.e. all outSubDirList[i] in outFolder
	for (i = 0; i < outSubDirList.length; i++) {
		if (File.exists(outFolder+outSubDirList[i]) == false){
			File.makeDirectory(outFolder+outSubDirList[i]);
		}
	}
	
	// load all stacks and proceed through manual tracing
	print("loading stacks for manual egg chamber annotation...");
	for (i = 0; i < fileList.length; i++) {
		//open current file in the list
		curFileName = inFolder+outSubDirList[i]+fileList[i];
		print(" ");
		print("opening file "+curFileName);
		open(curFileName);
		
		// let user circle the egg chambers manually
		roiManager("reset");
		waitForUser("Annotate stack manually", 
		"Circle each of the egg chambers in the stack using the free selection tool, \n"
		+ "then add each selection to the ROI manager (CMD+t).\n"
		+ "A nucleus should overlap >90% with the egg chamber ROI for correct assignment.\n"
		+ "It's ok if ROIs overlap a tiny bit. \n"
		+ "Click ok when you are done generating ROIs for all egg chambers in this stack.");
		
		// convert each egg chamber ROI into a 2D mask and
		// save in dedicated subfolder of the output dir.
		EggChamberSegFolderName = "eggChamberSEG/";
		saveEggChamberROIs(outFolder,outSubDirList[i],fileList[i],EggChamberSegFolderName);	
		
		// close stack
		close();		
	}	
	print("Done.");
}

function saveEggChamberROIs(outFolder,outSubDir,originalImgTitle,EggChamberSegFolderName) { 

	selectWindow(originalImgTitle);	
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	
	// collect the name of the image without the file extension to use as an output subfolder
	k = lastIndexOf(originalImgTitle, ".");
	if(k!=-1){
		imgNameWOExt = substring(originalImgTitle,0,k);
	}else{
		imgNameWOExt = originalImgTitle;
	}
		
	// build output subdirectory with the same name as the image (if needed)
	eggChamberDir = outFolder+outSubDir+imgNameWOExt+"/";
	if (File.exists(eggChamberDir) == false){
			File.makeDirectory(eggChamberDir);
	}
	
	// build eggChamberTif subdirectory inside the output subdirectory just created (if needed)
	saveDir = eggChamberDir+EggChamberSegFolderName;
	if (File.exists(saveDir) == false){
			File.makeDirectory(saveDir);
	}
	
	// loop through ROIs in ROI manager and generate a 2D mask image for each ROI
	// where the values inside the mask encode the index of the ROI.
	// (adding offset of 1 so that mask values start at 1, not 0).
	roiManager("show all");
	n = roiManager("count");
	print("Saving "+n+" egg chamber masks...");
	for (i = 0; i < n; i++) {
		selectWindow(originalImgTitle);
		roiManager("Select", i);
		run("Create Mask");
		run("16-bit");
		run("Divide...", "value=255");
		idx = i+1; // adding offset so that output ROI masks start at 1, not 0
		run("Multiply...", "value="+idx);
		setMinAndMax(0, idx);
		
		// save and close egg chamber segmentation mask
		imgName = imgNameWOExt+"_eggChamber"+idx+".tif";
		rename(imgName);
		save(saveDir+imgName);
		close();
	}	
	return;
}
