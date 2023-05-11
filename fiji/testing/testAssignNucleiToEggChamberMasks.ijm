macro "test assignNucleiToEggChamberMasks"{
	
	imgNameWOExt = "8-Ctrl-646,MPM2-488,Ser5ph-Cy3-1zCorr";
	resultImgName = "res";
	outFolder = "/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/";
	outSubDir = "Ctrl/";
	EggChamberSegFolderName = "eggChamberSEG/";
	
	// open image stack that includes nuclei as last channel
	open("/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/Ctrl/"
	+"8-Ctrl-646,MPM2-488,Ser5ph-Cy3-1zCorr/eggChamberTIF/"
	+"8-Ctrl-646,MPM2-488,Ser5ph-Cy3-1zCorrZcorrFinalNucMask.tif");
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	run("Duplicate...", "title=tmpNucMasks duplicate channels=C");
	inputNucMasks = getTitle();
	assignNucleiToEggChamberMasks(inputNucMasks,imgNameWOExt,resultImgName,
											outFolder,outSubDir,EggChamberSegFolderName);
}


// Function which finds which egg chamber each nucleus belongs to and
// stores that info into a z-stack of the nuclei masks, where the value in each
// nucleus mask is the egg chamber ID of the nucleus.
// (this is not efficient at all in terms of data storage - a hash table or 
// similar dictionary would be sufficient, but the stack format makes it very easy down the 
// road during the analysis steps.)

// Load egg chamber stack with nuclei segmentation (or already loaded)
// Locate matching eggchamberSEG folder x
// Get list of eggChamber mask files x

// Build empty stack that will store the eggChamber IDs. x
// Loop through segmented nuclei (find their number) 
//	Copy nucleus mask stack, isolating only the current nucleus and set its value to one
//	Max project nucleus 1-mask 
//	Compute sum of voxels in image (ie. Nucleus surface in pixel units)
//	Compute sum of voxels in image (ie. Nucleus surface in pixel units)
//  set eggChamberFound flag to 0
//  	Loop through egg chamber mask files x
//  	Load 2D mask. x
//		Multiply by eggchamber mask
//		Verify whether sum of voxel of multiplied image is greater than <0.9>*<nuc surface>*<eggChamber ID>
//		build 2-column array listing nuclei IDs and corresponding egg chamber IDs.
//		if nucleus overlaps >90% with eggchamber
//			multiply single nucleus 1-mask stack by egg chamber ID.
//			add this stack to the stack storing all nuclei egg chamber IDs.
//			set eggChamberFound flag to 1, to move to next nucleus

//output egg chamber nucleus IDs channel to data stack; 
function assignNucleiToEggChamberMasks(inputNucMasks,imgNameWOExt,resultImgName,
											outFolder,outSubDir,EggChamberSegFolderName){
	
	// sets the minimum fraction of the nucleus mask (its max projection) that has to overlap with the
	// egg chamber 2D mask in order to assign the nucleus. needs to be between 0 and 1, recommended value 0.9
	minOverlapThreshold = 0.8;
	
	// check whether the egg chamber segmentation folder exists.
	eggChamberDir = outFolder+outSubDir+imgNameWOExt+"/"+EggChamberSegFolderName;
	if (File.exists(eggChamberDir) == false){
		print("Missing egg Chamber segmentation folder: "+ eggChamberDir);
		print("Skipping...");
		return;
	}
	
	// list image files within the egg chamber segmentation folder 
	// collect names of image files in the input dir and its subfolders and
	// store into outSubDirList & fileList, so that the path to file i is
	// inFolder+outSubDirList[i]+fileList[i].
	dirList = getFileList(eggChamberDir);
	fileList = newArray(1000);
	outSubDirList = newArray(1000);
	ctr = 0;
	print(" ");
	for (i = 0; i < lengthOf(dirList); i++) {
	    if (endsWith(dirList[i], ".tif")) { 
	        fileList[ctr] = dirList[i];
	        ctr = ctr+1;
	    } 
	}
	fileList = Array.trim(fileList, ctr);
	
	if(ctr == 0){
		print("No files found in egg Chamber segmentation folder: "+ eggChamberDir);
		print("Skipping...");
		return;
	}
	
	// build a stack that will replicate the nuclei masks, 
	// but where each nucleus value encodes the ID of the egg chamber assigned. 
	selectWindow(inputNucMasks);	
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	newImage("eggChamberIDs", "16-bit", sizeX, sizeY, sizeZ);
	
	// find number of nuclei in z-stack
	selectWindow(inputNucMasks);
	run("Z Project...", "projection=[Max Intensity]");
	run("Select All");
	getStatistics(area, mean, min, max, std);
	maxNucleiFound = max;
	close(); // close max projection
	
	// make sure that integrated Density is part of the measurements
	run("Set Measurements...", "area mean standard min" +
	" integrated median area_fraction display redirect=None decimal=3");
	
	// loop through all putative nuclei
	for(idx = 1; idx <= maxNucleiFound; idx++){
		
		// check whether nucleus with ID = idx exists
		selectWindow(inputNucMasks);
		run("Duplicate...", "duplicate");
		// create stack duplicateIdx which will be 0 everywhere, 
		// except where the current nucleus (idx) is which is set to 1.
		rename("duplicateIdx"); 
		setAutoThreshold("Default dark");
		setThreshold(idx, idx, "raw");
		run("Convert to Mask", "method=Default background=Dark black");
		run("16-bit");
		run("Divide...", "value=255.000 stack");
		run("Z Project...", "projection=[Max Intensity]");
		rename("maxDuplicateIdx");
		getStatistics(area, mean, min, max, std, histogram);
		if(max == 1){
			nucleusFound = 1;
		}else {
			nucleusFound = 0;
		}
		
		// if nucleus exists, compute the max egg chamber overlap on the max projection
		if(nucleusFound == 1){
			
			// loop through all egg chamber 2D masks present in the folder
			maxEggChamberOverlap = 0;
			maxEggChamberID = 0;
			
			// compute surface area of nucleus max projection
			selectWindow("maxDuplicateIdx");
			run("Select All");
			run("Measure");
			nucArea = getResult("IntDen", nResults-1);
			
			for(i = 0; i < fileList.length; i++){
				
				curFileName = eggChamberDir+fileList[i];
				open(curFileName);	
				rename("curEC");
				getStatistics(area, mean, min, max, std, histogram);
				curECval = max;
				print("curECval = "+curECval);
				// check that the loaded img has expected dimensions
				getDimensions(ecSizeX, ecSizeY, ecC, ecSizeZ, F);
				
				// compute overlap
				if ((ecSizeX == sizeX) && (ecSizeY == sizeY) && (ecSizeZ == 1)){
					imageCalculator("Multiply create", "curEC","maxDuplicateIdx");
					rename("tmpMult");
					run("Select All");
					run("Measure");
    				curOverlap = getResult("IntDen", nResults-1);
    				print("curOverlap = "+curOverlap);
    				m = curOverlap/(curECval*nucArea);
    				print("idx = "+idx+"i = "+i+"; curECval = "+curECval+"; curOverlap = "+curOverlap
    					+"; m = "+m+"; maxEggChamberOverlap = "+maxEggChamberOverlap);
    				if(m > maxEggChamberOverlap){
    					maxEggChamberOverlap = m;
    					maxEggChamberID = curECval;
    				}
    				close("tmpMult");
				}
				close("curEC");
			}
			
			// multiply current nucleus 3D mask by egg chamber value that maximized the overlap
			// and add to place holder z-stack
			selectWindow("duplicateIdx");
			if(maxEggChamberOverlap > minOverlapThreshold){
				print("maxEggChamberID = "+maxEggChamberID);
				run("Multiply...", "value="+maxEggChamberID+" stack");
				imageCalculator("Add stack", "eggChamberIDs","duplicateIdx");
			}
		}
		close("duplicateIdx");
		close("maxDuplicateIdx");
	}
	
	// insert the egg chamber segmentation masks as a channel right before the channel where the nuclei masks were
	selectWindow("eggChamberIDs");	
	rename(resultImgName);
}

// appends an extra color channel (channelSource) to a hyperstack (imgSource) and renames the result newImgName 
function addChannelToImg(imgSource,channelSource,newImgName,keepSourceImgs){

	selectWindow(imgSource);
	b1 = bitDepth();
	selectWindow(channelSource);
	b2 = bitDepth();
	if ((b1 == 32) | (b2 == 32)){
		selectWindow(imgSource);
		run("32-bit");
		selectWindow(channelSource);
		run("32-bit");
	}
	
	selectWindow(imgSource);
	b1 = bitDepth();
	if(keepSourceImgs == 1){
		keepString = " keep";
	}else{
		keepString = "";
	}
	
	selectWindow(imgSource);
	getDimensions(w, h, c, z, f);
	if(c==1){
		argumentString1 = "c1="+imgSource+" c2="+channelSource+" create"+keepString;
		argumentString0 = "c1="+imgSource+" c2="+channelSource+" create";	
	}else{
		run("Split Channels");
		argumentString0 = "";
		for(i=1;i<=c;i++){
			argumentString0 = argumentString0 + "c"+i+"=C"+i+"-"+imgSource+" ";
		}
		i = c+1;
		argumentString1 = argumentString0 + "c"+i+"="+channelSource+" create"+keepString;
		argumentString0 = argumentString0 + " create";
	}
	
	if(keepSourceImgs == 1){
		// add channel to split channels from original image, keep sources
		run("Merge Channels...", argumentString1+keepString);
		rename(newImgName);

		if(c!=1){
			//recombine split channels into original image
			run("Merge Channels...", argumentString0);
			rename(imgSource);
		}
	}else{
		run("Merge Channels...", argumentString1);
		rename(newImgName);
	}
}