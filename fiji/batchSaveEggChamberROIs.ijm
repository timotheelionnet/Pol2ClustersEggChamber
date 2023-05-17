// this macro runs through a folder of raw data and loads each stack sequentially,
// prompting the user to segment egg chambers manually. This is done using the free
// drawing tool to generate ROIs and then enter them for each stack into the ROI manager.
// the macro then saves the ROIs as 2D masks in dedicated subfolders for each image.

// tracing tips: use the plane where the egg chamber is the widest and trace the boundaries. 
// it's ok if egg chambers overlap a bit. The idea is that a nuclei will be assigned to an
// egg chamber if > 80 % of its volume falls within the 2D mask (extended vertically).

// expected input folder structure:
// <inFolder/> 
//		|__ <condition1/>
//				<sample1.tif>	
//				<sample2.tif>					
//				...
//		|__ <condition2/>
//				<sample1.tif>	
//				<sample2.tif>				
//				...
// the script should also process images located directly in inFolder (i.e. not in a <condition> subfolder),
// but this hasn't been tested. SubSubFolders are NOT supported.

// output folder structure corresponding to the above input folder structure:
// <outFolder/> 
//		|__ <condition1/>
//				|__<sample1/>
//						|__ <eggChamberSEG/>	
//								|__sample1_eggChamberAnnotations.csv
//								|__sample1_eggChamber1.tif	
//								|__sample1_eggChamber2.tif	
//								...	
//				|__<sample2/>
//						|__ <eggChamberSEG/>	
//								|__sample2_eggChamberAnnotations.csv
//								|__sample2_eggChamber1.tif	
//								|__sample2_eggChamber2.tif	
//								...	
//				|__ ...
//
//		|__ <condition2/>
//				|__<sample1/>
//						|__ <eggChamberSEG/>	
//								|__sample1_eggChamberAnnotations.csv
//								|__sample1_eggChamber1.tif	
//								|__sample1_eggChamber2.tif	
//								...	
// 		|__ ...

// each output tif file is a 2D mask of one of the egg chamber masks and the value
// inside the mask matches the index of the egg chamber in the file name.
// the "Annotations.csv" file for each egg chamber is a place holder csv to use
// as a means to later manually annotate the development stages of each egg chamber. 
// 2 columns file, first column is the index of the egg chamber, second column is initialized to zeros.
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
	print("Starting batchSaveEggChamberROIs ...");
	for (i = 0; i < lengthOf(dirList); i++) {
	    if (endsWith(dirList[i], ".tif") || endsWith(dirList[i], ".nd2")) { 
	        fileList[ctr] = dirList[i];
	        outSubDirList[ctr] = "";
	        print("fname: "+fileList[ctr]+"; dir: "+outSubDirList[ctr]+"; ctr: "+ctr);
	        ctr = ctr+1;
	    } else{
	    	subDirList = getFileList(inFolder+dirList[i]);
	    	if (subDirList.length>0){
	    		for (j = 0; j < lengthOf(subDirList); j++) {
	    			if (endsWith(subDirList[j], ".tif") || endsWith(subDirList[j], ".nd2")) { 
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
	
	// create analysis subfolders for each experimental condition (if not already there), 
	// i.e. all outSubDirList[i] in outFolder
	for (i = 0; i < outSubDirList.length; i++) {
		if (File.exists(outFolder+outSubDirList[i]) == false){
			File.makeDirectory(outFolder+outSubDirList[i]);
		}
	}
	
	// create analysis sub-subfolders for each image (if not already there)
	for (i = 0; i < outSubDirList.length; i++) {
		k = lastIndexOf(fileList[i], ".");
		if(k!=-1){
		    imgNameWOExt = substring(fileList[i],0,k);
		}else{
		    imgNameWOExt = fileList[i];
		}
		eggChamberDir = outFolder+outSubDirList[i]+imgNameWOExt+"/";
		if (File.exists(eggChamberDir) == false){
		        File.makeDirectory(eggChamberDir);
		}
	}
	
	// create eggChamberSEG sub-subfolders in each analysis subfolder if needed
	// give warning if it already exists.
	eggChamberSegFolderName = "eggChamberSEG/";
	warningCleared = 0;
	for (i = 0; i < outSubDirList.length; i++) {
		k = lastIndexOf(fileList[i], ".");
		if(k!=-1){
		    imgNameWOExt = substring(fileList[i],0,k);
		}else{
		    imgNameWOExt = fileList[i];
		}
		eggChamberDir = outFolder+outSubDirList[i]+imgNameWOExt+"/";
		saveDir = eggChamberDir+eggChamberSegFolderName;
		if (File.exists(saveDir) == false){
		        File.makeDirectory(saveDir);
		}else{
			if(warningCleared == 0){
				waitForUser("Analysis Folder already exists!", 
					"Warning: this analysis will erase ALL prior egg chamber masks for the entire dataset\n"
					+"Click OK to proceed, CANCEL ortherwise");
				warningCleared = 1;
			}
			deleteList = getFileList(saveDir);
			for(j=0;j<deleteList.length;j++){
				File.delete(saveDir+deleteList[j]);
			}
			print("deleted " + saveDir);
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
		+ "A nucleus should overlap >80% with the egg chamber ROI for correct assignment.\n"
		+ "It's ok if ROIs overlap a tiny bit. \n"
		+ "DO NOT HIT OK now - do it only once you are done generating ROIs for all egg chambers in this stack!");

		if(roiManager("count") == 0){
			waitForUser("No Annotation Found!", "It looks like you did not add any egg chambers to the ROI selector!\n" 
				+"Click OK if you want to proceed without saving any egg chamber from this stack.\n"
				+"Otherwise circle each of the egg chambers using the free selection tool, \n"
				+ "then add each selection to the ROI manager (CMD+t).\n"
				+"Click OK when finished to save.");
		}
		// convert each egg chamber ROI into a 2D mask and
		// save in dedicated subfolder of the output dir.
		saveEggChamberROIs(outFolder,outSubDirList[i],fileList[i],eggChamberSegFolderName);	
		
		// close stack
		close();	
		
		// generate in the output folder a CSV file that will hold in a first column
		// all the egg chamber indices and for each its corresponding developmental stage
		// in a second column
		// (developmental stages left blank at this stage, use the file as a template)
		eggChamberSuffix = "_eggChamber";
		generateEggChamberCSV(outFolder,outSubDirList[i],
					fileList[i],eggChamberSegFolderName,eggChamberSuffix);
	}	
	print("Done.");
}

// generate the name of the analysis subfolder with the name of the original z-stack
function generateZStackOutputSubfolderName(outFolder,outSubDir,originalImgTitle){
	// collect the name of the image without the file extension to use as an output subfolder
	k = lastIndexOf(originalImgTitle, ".");
	if(k!=-1){
		imgNameWOExt = substring(originalImgTitle,0,k);
	}else{
		imgNameWOExt = originalImgTitle;
	}
		
	// build output subdirectory with the same name as the image (if needed)
	eggChamberDir = outFolder+outSubDir+imgNameWOExt+"/";
	
	return eggChamberDir;
}

// generate the name of the subfolder of the analysis subfolder called eggChamberSEG
function generateEggChamberSegOutputSubfolderName(outFolder,outSubDir,originalImgTitle,EggChamberSegFolderName){
	// collect the name of the image without the file extension to use as an output subfolder
	k = lastIndexOf(originalImgTitle, ".");
	if(k!=-1){
		imgNameWOExt = substring(originalImgTitle,0,k);
	}else{
		imgNameWOExt = originalImgTitle;
	}
		
	// build output subdirectory with the same name as the image (if needed)
	eggChamberDir = outFolder+outSubDir+imgNameWOExt+"/";
	
	// build eggChamberTif subdirectory inside the output subdirectory just created (if needed)
	saveDir = eggChamberDir+EggChamberSegFolderName;
	
	return saveDir;
}


// looks through the output folder for egg chamber segmentations
// and generates one csv file that lists each egg chamber index (col 1)
// and the a placeholder second column destined to hold the developmental stages (col 2)
// to be manually curated later.
function generateEggChamberCSV(outFolder,subDirName,
						fileName,eggChamberSegFolderName,eggChamberSuffix){
	
	// saving directory name
	ecDir = generateEggChamberSegOutputSubfolderName(
		outFolder,subDirName,fileName,eggChamberSegFolderName);
	
	// find all the egg chamber files in the folder and extract their indices into the 
	// array ecIdx
	fList =  getFileList(ecDir);
	ctr = 0;
	ecIdx = newArray(lengthOf(fList));
	for(i=0;i<lengthOf(fList);i++){
		if(fList[i].contains(eggChamberSuffix)){
			extIdx = lastIndexOf(fList[i], ".");
			ecIdx[ctr] = substring(fList[i], 
						indexOf(fList[i],eggChamberSuffix)+lengthOf(eggChamberSuffix), 
						extIdx);			
			ctr = ctr+1;			
		}
	}
	ecIdx = Array.trim(ecIdx,ctr);
	
	// generate a text file in the eggChamber dir that will store the indices of
	// each egg chamber and a place holder second column for the stage
	csvFilePath = ecDir+"eggChamberStages.csv";
	placeHolderColumn = newArray(ctr);
	generate2columnCSVfromArrays(csvFilePath,ecIdx,placeHolderColumn,
								"eggChamberID","eggChamberStage");
	
}

// saves 2 arrays as a csv file in 2 columns. Arrays need to be the same size.
function generate2columnCSVfromArrays(filePath,col1Array,col2Array,header1,header2){
    
    if (lengthOf(col1Array) != lengthOf(col2Array)){
    	print("2 arrays have different sizes, cannot save as csv table");
    	return;
    }
    // Create an empty string to store the CSV data
    csvData = header1+","+header2+"\n";
    // Loop through the numbers 1 to n and add them to the CSV data with a zero in the second column
    for (i = 0; i < lengthOf(col1Array); i++) {
        csvData = csvData + col1Array[i] + ","+col2Array[i] + "\n";
    }
    
    // Save the CSV file
   f = File.open(filePath);
   print(f,csvData);
   File.close(f);
}

// Loops through the egg chamber outline ROIs stored in the ROI manager and for each: 
// 	- generates a 2D mask image with the value corresponding to the current ROI ID (ROI 0: value 1; ROI 1: value 2, etc)
// 	- saves the 2D mask images in a subfolder called EggChamberSegFolderName

function saveEggChamberROIs(outFolder,outSubDir,originalImgTitle,eggChamberSegFolderName) { 

	selectWindow(originalImgTitle);	
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	
	// get the name of the output eggChamberTif subdirectory 
	saveDir = generateEggChamberSegOutputSubfolderName(
		outFolder,outSubDir,originalImgTitle,eggChamberSegFolderName);
	
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
		run("Grays");
		
		// save and close egg chamber segmentation mask
		imgName = imgNameWOExt+"_eggChamber"+idx+".tif";
		rename(imgName);
		save(saveDir+imgName);
		close();
	}	
	return;
}
