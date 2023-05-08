var defaultHoechstChannel = 1;
var defaultAvgNucleiDiameterInUm = 20;
var defaultSaveInitialSegResults = true;
var defaultUseMorphologyFilters = true;
var defaultMinVolume =  400;
var defaultMaxVolume =  40000;
var defaultMaxSurfToVolRatio = 1;
var defaultMinSphericity =  0.4;
var defaultMaxKurtosis = 3;
var defaultMinCV =  0.15;
var defaultMaxCV =  0.45;

var hoechstChannel;
var avgNucleiDiameterInUm;
var saveInitialSegResults;
var useMorphologyFilters;
var minVolume;
var maxVolume;
var maxSurfToVolRatio;
var minSphericity;
var maxKurtosis;
var minCV;
var maxCV;
var useEggChamberMasks;

macro "segmentNuclei"{
	run("Close All");
	setBatchMode(true);
	inputParameters();

	inFolder = getDirectory("choose the input directory");
	outFolder = getDirectory("choose the output directory");
	
	//inFolder = "/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/test/";
	//outFolder = "/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/junkOut/";
	
	// collect names of image files in input dir and its subfolders and
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
	
	print("segmenting nuclei...");
	for (i = 0; i < fileList.length; i++) {
		//open current file in the list
		curFileName = inFolder+outSubDirList[i]+fileList[i];
		print(" ");
		print("opening file "+curFileName);
		open(curFileName);
		
		//remove extension from filename
		originalImgTitle = getTitle();	
		if(indexOf(originalImgTitle,'.')>0){
			imgNameWOExt = originalImgTitle.substring(0,indexOf(originalImgTitle,'.'));
		}else{
			imgNameWOExt = originalImgTitle;
		}
		rename(imgNameWOExt);
		originalImgTitle = imgNameWOExt;
		
		// de-trend intensity variation along z - output is an image called zcorrImgTitle 
		// which contains the corrected data
		print("de-trending intensity along Z...");
		zcorrImgTitle = originalImgTitle+"zCorr";
		eggChamberIntensityMeasurementAllChannels2("eggChamber",
			originalImgTitle,hoechstChannel,zcorrImgTitle);
		close(originalImgTitle);
		run("Collect Garbage");
			
		// segment nuclei - output is a single channel z-stack
		// called initNucMasksTitle which contains the nuclei masks
		print("initial nuclei segmentation...");
		initNucMasksTitle = "initNucleiMasks";
		segmentNuclei3D(zcorrImgTitle,initNucMasksTitle,hoechstChannel,avgNucleiDiameterInUm);
		run("Collect Garbage");
		
		// save initial segmentation if option was selected.
		if(saveInitialSegResults == true){
			print("saving initial nuclei segmentation...");
			addChannelToImg(zcorrImgTitle,initNucMasksTitle,"initNucOut",1);
			EggChamberTifFolderName = "eggChamberTIF/";
			fileSuffix = "zCorrInitNucMask.tif";
			saveInitSegResult("initNucOut",outFolder,outSubDirList[i],
				EggChamberTifFolderName,fileList[i],fileSuffix);
			print("done saving");
			close("initNucOut");
			run("Collect Garbage");
		}
		
		//clean up aberrantly segmented objects based on morphological filters
		// output is an image called finalNucMasksTitle containing the updated masks
		print("cleaning up segmentation with morphological filters...");
		finalNucMasksTitle = "finalNucMasks";
		selectCorrectlySegmentedNuclei2(zcorrImgTitle,initNucMasksTitle,hoechstChannel,finalNucMasksTitle,
			minVolume,maxVolume, maxSurfToVolRatio, minSphericity,
			maxKurtosis,minCV,maxCV);
		close(initNucMasksTitle);
		run("Collect Garbage");
		
		// generate mask z-stack encoding the egg chamber assignment of each nucleus
		if(useEggChamberMasks == true){
			print("assigning nuclei to respective egg chambers...");
			EggChamberSegFolderName = "eggChamberSEG/";
			eggChamberIDsTitle = "eggChamberIDs";
			
			eggChamberDataFound = assignNucleiToEggChamberMasks(
											finalNucMasksTitle,imgNameWOExt,eggChamberIDsTitle,
											outFolder,outSubDirList[i],EggChamberSegFolderName);
			if(eggChamberDataFound == 1){
				// append to detrended image so the eggchamber ID is measured with the metrics								
				addChannelToImg(zcorrImgTitle,eggChamberIDsTitle,zcorrImgTitle,0);			
			}	
		}
		
		// compute and save metrics
		print("computing geometry and intensity metrics on objects...");
		EggChamberCsvFolderName = "eggChamberCSV/";
		runMetricsAndSave(finalNucMasksTitle,zcorrImgTitle,EggChamberCsvFolderName,
			outFolder,outSubDirList[i],fileList[i],"Geom.csv","Int.csv");
		
		// compute average intensities for each channel across entire image	for background subtraction 
		// during later data analysis steps	
		print("computing whole image intensity metrics...");	
		wholeImgIntensityMeasurementAllChannels(outFolder+outSubDirList[i]+ 
			originalImgTitle+"/"+EggChamberCsvFolderName,zcorrImgTitle);
		
		// compute average intensities for each channel across broader egg chambers region for background subtraction
		// during later data analysis steps		
		print("computing egg chamber intensity metrics...");	
		eggChamberIntensityMeasurementAllChannels(outFolder+outSubDirList[i]+ 
			originalImgTitle+"/"+EggChamberCsvFolderName,zcorrImgTitle,hoechstChannel);
		run("Collect Garbage");
		
		// merge seg result w de-trended z-stack and save
		print("saving segmentation results...");
		addChannelToImg(zcorrImgTitle,finalNucMasksTitle,"finalNucOut",1);
		saveFileName = imgNameWOExt+"ZcorrFinalNucMask.tif";
		selectWindow("finalNucOut");
		save(outFolder+outSubDirList[i]+imgNameWOExt+"/"+EggChamberTifFolderName 
			+ saveFileName);
		run("Close All");
		run("Collect Garbage");
	}
	setBatchMode("exit and display");
  	print("done.");
}

// Function which finds which egg chamber each nucleus belongs to and
// stores that info into a z-stack of the nuclei masks, where the value in each
// nucleus mask is the egg chamber ID of the nucleus.
// (this is not efficient in terms of data storage - a hash table or 
// similar dictionary would be better from a space perspective, 
// but the stack format makes it very easy down the road during the analysis steps.)

// input: nuclei segmentation stack (inputNucMasks)
// egg chamber segmentation 2D masks are expected to be located at the following path:
// 		outFolder/outSubDir/imgNameWOExt/EggChamberSegFolderName/

// resultImgName: image stack output as the result of the analysis 
// (original stack left intact at the end). 

// output is a 0/1 flag

function assignNucleiToEggChamberMasks(inputNucMasks,imgNameWOExt,resultImgName,
											outFolder,outSubDir,EggChamberSegFolderName){
	//initialize the output flag
	eggChamberDataFound = 1; 
	
	// this variable sets the minimum fraction of the nucleus mask (its max projection) 
	// that has to overlap with the egg chamber 2D mask in order to assign the nucleus. 
	// needs to be between 0 and 1, recommended value 0.9
	minOverlapThreshold = 0.9;
	
	// check whether the egg chamber segmentation folder exists.
	eggChamberDir = outFolder+outSubDir+imgNameWOExt+"/"+EggChamberSegFolderName;
	if (File.exists(eggChamberDir) == false){
		print("Missing egg Chamber segmentation folder: "+ eggChamberDir);
		print("Skipping...");
		eggChamberDataFound = 0;
		return eggChamberDataFound;
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
		eggChamberDataFound = 0;
		return eggChamberDataFound;
	}
	
	// build a stack that will replicate the nuclei masks, 
	// but where each nucleus value encodes the ID of the egg chamber assigned. 
	selectWindow(inputNucMasks);	
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	newImage("eggChamberIDs", "16-bit", sizeX, sizeY, sizeZ);
	
	// find number of nuclei in z-stack
	selectWindow(inputNucMasks);
	run("Z Project...", "projection=[Max Intensity]");
	getStatistics(area, mean, min, max, std);
	maxNucleiFound = max;
	close(); // close max projection
	
	// make sure that integrated Density is part of the measurements
	run("Set Measurements...", "area mean standard min" +
	" integrated median area_fraction display redirect=None decimal=3");
	
	// loop through all putative nuclei
	for(idx = 1; idx < maxNucleiFound; idx++){
		
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
	if(isOpen("Results")){
		close("Results");
	}
	return eggChamberDataFound;
}

// takes an image (inputWindowName) generates a global mask for the entire egg chamber
// using the intensity in measurementChannel as a basis (and using some filtering/thresholding), 
// then measures intensity/neighbors etc in all channels using the new mask as an ROI
// and saves the resulting tables in the folder savePath with names C1_eggChamberInt.csv, C2_eggChamberInt.csv etc
function eggChamberIntensityMeasurementAllChannels(savePath,inputWindowName,measurementChannel){
	// rescale by a factor of 4 for faster processing
	selectWindow(inputWindowName);
	run("Duplicate...", "duplicate");
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	x2 = round(sizeX/4);
	y2 = round(sizeY/4);
	z2 = round(sizeZ/4);
	run("Size...", "width="+x2+" height="+y2+" depth="+z2
		+" constrain average interpolation=Bicubic");
	rename("resized");
	
	// generate a mask that encompasses the egg chamber, to be used for measurements
	selectWindow("resized");
	run("Duplicate...", "duplicate channels="+measurementChannel);
	rename("tmpDuplicate1");
	selectWindow("tmpDuplicate1");
	run("Median...", "radius=4 stack");
	run("Auto Threshold", "method=Default white stack");
	setAutoThreshold("Default");
	run("Threshold...");
	setThreshold(255, 255);
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Light black");
	
	selectWindow(inputWindowName);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	for (i = 1; i <= C; i++) {
	    selectWindow("resized");
	    run("Duplicate...", "duplicate channels="+i);
		rename("tmpDuplicate2");
		
		run("Intensity Measurements 2D/3D", "input=tmpDuplicate2 labels=tmpDuplicate1"+
			" mean stddev max min median mode skewness kurtosis numberofvoxels volume");
			
		saveAs("tmpDuplicate2-intensity-measurements",savePath+"C"+i+"_eggChamberInt.csv");
		run("Close");
		
		selectWindow("tmpDuplicate2");
		close();
	}

	print("intensity over mask file saved");
	selectWindow("tmpDuplicate1");
	close();
	selectWindow("resized");
	close();
}

// appends an extra color channel (channelSource) to a hyperstack (imgSource) and renames the result newImgName 
function addChannelToImg(imgSource,channelSource,newImgName,keepSourceImgs){
	if(imgSource == newImgName){
		selectWindow(imgSource);
		rename("tmpImg1AddChannelToImg");
		imgSource = "tmpImg1AddChannelToImg";
	}
	if(channelSource == newImgName){
		selectWindow(channelSource);
		rename("tmpImg2AddChannelToImg");
		channelSource = "tmpImg2AddChannelToImg";
	}
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
		if(isOpen(channelSource)){
			close(channelSource);
		}
		if(isOpen(imgSource)){
			close(imgSource);
		}
	}
}

// performs measurements of intensity on the whole image
// savePath should just be a directory
function wholeImgIntensityMeasurementAllChannels(savePath,inputWindowName){
	selectWindow(inputWindowName);
	run("Duplicate...", "duplicate");
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	x2 = round(sizeX/4);
	y2 = round(sizeY/4);
	z2 = round(sizeZ/4);
	run("Size...", "width="+x2+" height="+y2+" depth="+z2
		+" constrain average interpolation=Bicubic");
	rename("resized");
	
	// generate a mask that encompasses the entire image, to be used for measurements
	run("Duplicate...", "duplicate channels=1");
	rename("tmpDuplicate1");
	selectWindow("tmpDuplicate1");
	setThreshold(0, 65535); 
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Light black");

	selectWindow(inputWindowName);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	for (i = 1; i <= C; i++) {
	    selectWindow("resized");
	    run("Duplicate...", "duplicate channels="+i);
		rename("tmpDuplicate2");
		run("Intensity Measurements 2D/3D", "input=tmpDuplicate2 labels=tmpDuplicate1"+
			" mean stddev max min median mode skewness kurtosis numberofvoxels volume");
		saveAs("tmpDuplicate2-intensity-measurements",savePath+"C"+i+"_wholeImgInt.csv");
		run("Close");
		
		selectWindow("tmpDuplicate2");
		close();
	}
	selectWindow("tmpDuplicate1");
	close();
	selectWindow("resized");
	close();
}

function selectCorrectlySegmentedNuclei2(inputImgTitle,maskImgTitle,
	hoechstChannel,outputImgTitle,
	minVolume,maxVolume, maxSurfToVolRatio,minSphericity,
	maxKurtosis,minCV,maxCV){
	
	// remove masks that touch any lateral boundary
	selectWindow(maskImgTitle);
	run("Remove Border Labels", "left right top bottom");
	rename(maskImgTitle);

	// collect table of object geometric features
	run("Analyze Regions 3D", "volume surface_area mean_breadth sphericity euler_number" 
	+" bounding_box centroid surface_area_method=[Crofton (13 dirs.)] euler_connectivity=26");

	// rename the morpholibj output table so ImageJ recognizes it as a Results table
	// the extension of the filename is not transfered to the morpholibJ table name
	if(indexOf(maskImgTitle,'.')>0){
		maskNameWOExt = maskImgTitle.substring(0,indexOf(maskImgTitle,'.'));
	}else{
		maskNameWOExt = maskImgTitle;
	}
	Table.rename(maskNameWOExt+"-morpho", "Results");
	
	// collect list of nuclei label that pass the geometry criteria
	labelToKeepGeom = newArray(nResults);
	gCtr = 0;
	for (i=0; i<nResults; i++){
		keepCurrent = 1;
		v = getResult("Volume",i);
		if((v<minVolume) || (v>maxVolume)){
			keepCurrent = 0;
		}
		stv = getResult("SurfaceArea",i)/getResult("Volume",i);
		if(stv > maxSurfToVolRatio){
			keepCurrent = 0;
		}
		s = getResult("Sphericity",i);
		if(s < minSphericity){
			keepCurrent = 0;
		}
		
		if(keepCurrent == 1){
			labelToKeepGeom[gCtr] = getResultString("Label",i);
			gCtr = gCtr+1;
		}
	}
	labelToKeepGeom = Array.trim(labelToKeepGeom, gCtr);
	close("Results");
	
	// compute intensity metrics on Hoechst channel
	selectWindow(inputImgTitle);
	run("Duplicate...", "title=hoechst duplicate channels="+hoechstChannel);
	run("Intensity Measurements 2D/3D", "input=hoechst labels="+ maskImgTitle +" mean stddev"
		+" max min median mode skewness kurtosis numberofvoxels volume");
	Table.rename("hoechst-intensity-measurements", "Results");
	
	// collect list of nuclei label that pass the intensity-based criteria
	labelToKeepInt = newArray(nResults);
	iCtr = 0;
	for (i=0; i<nResults; i++){
		keepCurrent = 1;
		
		k = getResult("Kurtosis",i);
		if(k > maxKurtosis){
			keepCurrent = 0;
		}
		cv = getResult("StdDev",i)/getResult("Mean",i);
		if((cv<minCV) || (cv>maxCV)){
			keepCurrent = 0;
		}
		
		if(keepCurrent == 1){
			labelToKeepInt[iCtr] = getResultString("Label",i);
			iCtr = iCtr+1;
		}
	}
	labelToKeepInt = Array.trim(labelToKeepInt, iCtr);
	close("Results");
	close("hoechst");
	
	// keep labels that satisfy both criteria
	labelToKeep = newArray(maxOf(iCtr,gCtr) );
	labelString = "label(s)=";
	ctr = 0;
	for(i=0;i<iCtr;i++){
		for (j = 0; j < gCtr; j++) {
			if(labelToKeepInt[i] == labelToKeepGeom[j]){
				labelToKeep[ctr] = labelToKeepInt[i];
				if(ctr == 0){
					labelString = labelString+labelToKeep[0];
				}else {
					labelString = labelString+","+labelToKeep[ctr];
				}
				ctr = ctr+1;
			}
		}
	}
	
	// remove non-selected nuclei from mask
	selectWindow(maskImgTitle);
	if (labelString != "label(s)="){
		run("Select Label(s)", labelString);
		print("Selected "+ ctr +" objects with following labels: "+labelString);
	}else{
		print("Found no objects to exclude.");
	}
	rename(outputImgTitle);
	return;
}

//computes geometry metrics, adds a few columns:
//1) file name;
//2) distance to image center 
// and saves the table in the file specified by savePath
function runMetricsAndSave(inputMaskTitle,inputDataTitle,EggChamberCsvFolderName,
	outFolder,outSubDir,inputFileName,geomSuffix,intSuffix){
	
	// remove extension from filename if needed 
	if(indexOf(inputFileName,'.')>0){
		inputFileNameWOExt = inputFileName.substring(0,indexOf(inputFileName,'.'));
	}else{
		inputFileNameWOExt = inputFileName;
	}
						
	// build output subdirectory with the same name as the image (if needed)
	eggChamberDir = outFolder+outSubDir+inputFileNameWOExt+"/";
	if (File.exists(eggChamberDir) == false){
			File.makeDirectory(eggChamberDir);
	}
	
	// build eggChamberCSV subdirectory inside the output subdirectory just created (if needed)
	saveDir = eggChamberDir+EggChamberCsvFolderName;
	if (File.exists(saveDir) == false){
			File.makeDirectory(saveDir);
	}
	
	selectWindow(inputMaskTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	run("Analyze Regions 3D", "volume surface_area mean_breadth"
		+" sphericity euler_number bounding_box centroid"
		+" surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");

	// rename the morpholibj output table so ImageJ recognizes it as a Results table
	// the extension of the filename is not transfered to the morpholibJ table name
	if(lastIndexOf(inputMaskTitle,'.')>0){
		imgNameWOExt = inputMaskTitle.substring(0,lastIndexOf(inputMaskTitle,'.'));
	}else{
		imgNameWOExt = inputMaskTitle;
	}
	print(imgNameWOExt);
	
	// computing distances to stack center
	Table.rename(imgNameWOExt+"-morpho", "Results");
	print("computing metrics across "+nResults+" objects...");
	
	// compute distance of each object centroid to the stack center
	getVoxelSize(dX, dY, dZ, u);
	iX = dX*sizeX/2;
	iY = dY*sizeY/2;
	iZ = dZ*sizeZ/2;
	
	distArray = newArray(nResults);
	for (i=0; i<nResults; i++){
	 	 oX = getResult("Centroid.X",i);
	 	 oY = getResult("Centroid.Y",i);
	 	 oZ = getResult("Centroid.Z",i);
	 	 
	 	 distArray[i] = sqrt( (oX-iX)*(oX-iX) + (oY-iY)*(oY-iY) + (oZ-iZ)*(oZ-iZ) );
	}
	
	// add distance to image center as a column to table
	Table.setColumn("Dist", distArray);
	
	fileNameArray = newArray(nResults);
	for (i=0; i<nResults; i++){
		fileNameArray[i] = inputFileName;
	}
	Table.setColumn("InputFileName", fileNameArray);
	
	//save and close geometry results
	saveAs("Results", saveDir + "allNuc" + geomSuffix);
	close("Results");
	
	// run intensity measurements on all channels
	selectWindow(inputDataTitle);
	getDimensions(w, h, c, z, f);
	for (i = 1; i <= c; i++) {
		selectWindow(inputDataTitle);
		run("Duplicate...", "duplicate channels="+i);
		rename("curChannel");
		run("Intensity Measurements 2D/3D", "input=curChannel"
		+" labels="+inputMaskTitle+" mean"
		+" stddev max min median mode skewness kurtosis numberofvoxels"
		+" volume");
		
		Table.rename("curChannel-intensity-measurements", "Results");
		saveAs("Results", saveDir + "C"+i+"_allNuc" + intSuffix);
		close("Results");
		close("curChannel");
	}
}

// self explanatory
function saveInitSegResult(imgName,outFolder,outSubDir,
			EggChamberTifFolderName,fileName,fileSuffix){
	
	// remove extension from filename if needed 
	if(indexOf(fileName,'.')>0){
		imgNameWOExt = fileName.substring(0,indexOf(fileName,'.'));
	}else{
		imgNameWOExt = fileName;
	}
						
	// build output subdirectory with the same name as the image (if needed)
	eggChamberDir = outFolder+outSubDir+imgNameWOExt+"/";
	if (File.exists(eggChamberDir) == false){
			File.makeDirectory(eggChamberDir);
	}
	
	// build eggChamberTif subdirectory inside the output subdirectory just created (if needed)
	saveDir = eggChamberDir+EggChamberTifFolderName;
	if (File.exists(saveDir) == false){
			File.makeDirectory(saveDir);
	}
	
	// save
	selectWindow(imgName);
	save(saveDir+imgNameWOExt+fileSuffix);
	return;
}


// self explanatory
function setParametersToDefault(){
	
	hoechstChannel = defaultHoechstChannel;
	avgNucleiDiameterInUm = defaultAvgNucleiDiameterInUm;
	useMorphologyFilters = defaultUseMorphologyFilters;
	minVolume = defaultMinVolume;
	maxVolume = defaultMaxVolume;
	maxSurfToVolRatio = defaultMaxSurfToVolRatio;
	minSphericity = defaultMinSphericity;
	maxKurtosis = defaultMaxKurtosis;
	minCV = defaultMinCV;
	maxCV = defaultMaxCV;
}

// self explanatory
function inputParameters(){
	useDefaults = true;
	
	setParametersToDefault();
	
	Dialog.create("Segmentation Parameters");
	
	Dialog.addMessage("******************* Initial nucleus Segmentation Module");
	Dialog.addNumber("Hoechst color channel:", hoechstChannel);	
	Dialog.addNumber("Average Nucleus Diameter in um ("+defaultHoechstChannel+"):", avgNucleiDiameterInUm);	
	Dialog.addCheckbox("Save initial segmentation results", true);
	
	Dialog.addMessage("******************* Morphology Filters Module (default values)");
	Dialog.addCheckbox("Use Morphology Filters", true);
	Dialog.addCheckbox("Use Default morphology parameters (will override any entry)", true);
	Dialog.addNumber("Min. Nucleus Volume in um^3 ("+defaultMinVolume+"):", minVolume);
	Dialog.addNumber("Max. Nucleus Volume in um^3 ("+defaultMaxVolume+"):", maxVolume);
	Dialog.addNumber("Max. Surface to Volume ratio ("+defaultMaxSurfToVolRatio+"):", maxSurfToVolRatio);
	Dialog.addNumber("Mim. Sphericity ("+defaultMinSphericity+"):",minSphericity);
	Dialog.addNumber("Max. Hoechst distribution Kurtosis ("+defaultMaxKurtosis+"):",maxKurtosis);
	Dialog.addNumber("Min. Hochst distribution CV ("+defaultMinCV+"):",minCV);
	Dialog.addNumber("Max. Hochst distribution CV ("+defaultMaxCV+"):",maxCV);
	
	Dialog.addMessage("******************* Egg Chamber masks");
	Dialog.addCheckbox("Use Egg Chamber 2D masks (requires them already generated)", true);
	Dialog.show();
	
	
	hoechstChannel = Dialog.getNumber();
	avgNucleiDiameterInUm = Dialog.getNumber();
	saveInitialSegResults = Dialog.getCheckbox();
	
	useMorphologyFilters = Dialog.getCheckbox();
	useDefaults = Dialog.getCheckbox();
	minVolume =  Dialog.getNumber();
	maxVolume =  Dialog.getNumber();
	maxSurfToVolRatio =  Dialog.getNumber();
	minSphericity =  Dialog.getNumber();
	maxKurtosis =  Dialog.getNumber();
	minCV =  Dialog.getNumber();
	maxCV =  Dialog.getNumber();
	useEggChamberMasks = Dialog.getCheckbox();
	
	if (useDefaults == true){
		setParametersToDefault();
	}
	
	print("*******************************");
	print("Hoechst color channel: "+hoechstChannel);
	print("Average Nuclei Diameter: ",avgNucleiDiameterInUm);
	print("Save initial segmentation results: ",saveInitialSegResults);
	print("use Morphology Filters: ",useMorphologyFilters);
	print("Use Default Morphology Filters: ",useDefaults);
	print("Min. nucleus Volume (um^2)",minVolume);
	print("Max. nucleus Volume (um^2)",maxVolume);
	print("Max. Surface-to-Volume Ratio: ",maxSurfToVolRatio);
	print("Min. Sphericity: ",minSphericity);
	print("Max. Hoechst distribution Kurtosis: ",maxKurtosis);
	print("Min. Hoechst distribution Coeff of Variation",minCV);
	print("Max. Hoechst distribution Coeff of Variation",maxCV);
}



// takes an input hyperstack inputWindowName, then 
// 1) subtracts the global minimum intensity in all channels (~black level offset correction)
// 2) generates a mask of the eggChamber using the intensity 
// 	in measurementChannel as a threshold (suggest using DAPI/Hoechst)
// 2) computes the trend of the mean intensity (computed within the mask) vs Z in each channel
// 3) fits the trend as a function of Z to a line in each channel
// 4) divides each channel intensity by the linear trend
// 5) outputs the result as deTrendedIntTitle
// maskType flag can be "eggChamber" to use the mask based on measurementChannel, or anything else, 
// in which case the whole image is used to compute the z trend.
function eggChamberIntensityMeasurementAllChannels2(maskType,
	inputWindowName,measurementChannel,deTrendedIntTitle){

	// generate a mask that encompasses the egg chamber (or the whole image, 
	// depending on the flag maskType), to be used for measurements of the avg intensity vs z.
	selectWindow(inputWindowName);
	run("Duplicate...", "duplicate channels="+measurementChannel);
	if (maskType == "eggChamber"){
		run("Bandpass Filter...", "filter_large=2000 filter_small=100 suppress=None tolerance=5 process");
		run("Auto Threshold", "method=Default white stack");
	}else{
		setThreshold(0, 65535);
		run("Convert to Mask", "method=Default background=Dark black");
	}
	rename("mask");
	
	// compute minimum intensity in each channel across the entire z-stack 
	//- this offset will then be subtracted from the raw data.
	computeMinInt(inputWindowName,"minWholeImg");

	// subtract channel minimum (computed across entire stack) from each channel 
	subtractBackgroundGlobally(inputWindowName,"minWholeImg","dataMinusOffset");
	
	// compute mean intensity within the mask in each slice z,c 
	// output is "meanEggChamber", a 16-bit 2D image (nz nc) holdig the mean
	computeMeanZTrend("dataMinusOffset","mask","meanEggChamber");
	
	// fit the the eggchamber mean intensity as a function of Z to a line
	// in each channel using "meanEggChamber" as an input.
	// output the linear trend "linTrend", a 16-bit 2D image (nz nc) holdig the fit
	computeLinearTrendAlongX("meanEggChamber","linTrend");

	// from input hyperstack dataMinusOffset and the 2D image "linTrend" holding the linear trend vs z by channel 
	// generates a corrected hyperstack deTrendedIntTitle (nx ny nz nc) where the intensity in each slice is corrected
	// Icorr(x,y,z,c) = Iinput(x,y,z,c) * linTrend(0,c) / linTrend(z,c)
	correctIntensityForZTrend("dataMinusOffset","linTrend",deTrendedIntTitle);

	selectWindow("linTrend");
	close();
	selectWindow("minWholeImg");
	close();
	selectWindow("meanEggChamber");
	close();
	selectWindow("dataMinusOffset");
	close();
	selectWindow("mask");
	close();

	// concert back to 16 bit
	selectWindow(deTrendedIntTitle);
	for (i = 1; i <= 4; i++) {
		Stack.setPosition(i, 1, 1);
		setMinAndMax(0, 65535);
	}
	run("16-bit");
}

// computes min intensity across each z,c slice of the input hyperstack
// returns it as an (nz,nc) 2D image
function computeMinInt(inputWindowName,minTitle){
	selectWindow(inputWindowName);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	newImage(minTitle, "16-bit", sizeZ, C, 1);
	for(izs = 0; izs<sizeZ; izs++) {
		for (i = 0; i < C; i++) {
		    selectWindow(inputWindowName);
		    Stack.setPosition(i+1, izs+1, 1);
		    run("Select All");
		    getStatistics(area, mean, minImg, max, std, histogram);
			selectWindow(minTitle);
			setPixel(izs, i, minImg);
		}
	}
}

// from an input hyperstack inputImgTitle (nx ny nz nc), 
// and minTitle, an nz,nc 2D image holding the minimum value at each slice
// computes the global minimum across the entire stack for each channel 
// and subtracts it from the input hyperstack. 
// Outputs the subtracted image called outputImgTitle
function subtractBackgroundGlobally(inputImgTitle,minTitle,outputImgTitle){
	selectWindow(inputImgTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	//duplicate input hyperstack
	run("Select All");
	run("Duplicate...", "duplicate");
	rename(outputImgTitle);
	run("32-bit");

	// compute absolute minimum in each channel, i.e. find the minimum of each row
	// of minTitle
	minArray = newArray(sizeY);
	selectWindow(minTitle);
	for (j = 0; j < C; j++) {
		for (i = 0; i < sizeZ; i++) {
			if(i==0){
				minArray[j] = getPixel(i,j);
			}else{
				x = getPixel(i,j);
				if(x<minArray[j]){
					minArray[j] = x;
				}
			}
		}
	}
	
	// subtract absolute minimum to each slice
	for(izs = 0; izs<sizeZ; izs++) {
		for(ic = 0; ic<C; ic++) {
			selectWindow(outputImgTitle);
			Stack.setPosition(ic+1, izs+1, 1);
			run("Subtract...", "value="+ minArray[ic] +" slice");
		}
	}
}	

// from an input hyperstack inputImgTitle (nx ny nz nc), 
// and a 2D image linTrendTitle (nz,nc) holding the linear trend vs z by channel 
// generates a corrected hyperstack outputImgTitle (nx ny nz nc) where the intensity in each slice is corrected
// Icorr(x,y,z,c) = Iinput(x,y,z,c) * linTrend(0,c) / linTrend(z,c)
function correctIntensityForZTrend(inputImgTitle,linTrendTitle,outputImgTitle){
	selectWindow(inputImgTitle);
	run("32-bit");
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	selectWindow(linTrendTitle);
	run("32-bit");
	getDimensions(trendX, trendY, trendC, trendZ, trendF);

	// check that dimensions of the inputs are consistent with each other
	if(trendX != sizeZ){
		print("trend has "+trendX+" rows, whereas the number of slices is "+sizeZ);
		return;
	}
	if(trendY != C){
		print("trend has "+trendY+" rows, whereas the number of channels is "+C);
		return;
	}

	// get the reference level in ech channel at the first frame
	trendRef = newArray(C);
	for(ic = 0; ic<C; ic++) {	
		selectWindow(linTrendTitle);
		trendRef[ic] = getPixel(0, ic);
	}

	// normalize each channel to its first z-slize
	selectWindow(inputImgTitle);
	run("Select All");
	run("Duplicate...", "duplicate");
	rename(outputImgTitle);
	run("32-bit");
	curTrend = newArray(C);
	for(izs = 0; izs<sizeZ; izs++) {
		for(ic = 0; ic<C; ic++) {
			selectWindow(linTrendTitle);
			curTrend[ic] = getPixel(izs, ic)/trendRef[ic];
			//print("ic="+ic"; izs = "+izs+"; curTrend: "+curTrend[ic]);
			
			selectWindow(outputImgTitle);
			Stack.setPosition(ic+1, izs+1, 1);
			run("Divide...", "value="+ curTrend[ic] +" slice");
		}
	}
}

// from hyperstack hsTitle (nx ny nz nc), 
// and a mask maskTitle (nx nz nz 1 8-bit stack) at each z-slice
// computes the mean within the mask as a function of z and c
// outputs is a 16-bit 2D image meanTrendTitle (nz nc) holdig the mean
function computeMeanZTrend(hsTitle,maskTitle,meanTrendTitle){
	selectWindow(hsTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	newImage(meanTrendTitle, "16-bit", sizeZ, C, 1);
	
	for(izs = 0; izs<sizeZ; izs++) {
		roiManager("reset");
	    selectWindow(maskTitle);
	    Stack.setPosition(1, izs+1, 1);
	    setThreshold(100, 255);
	    run("Create Selection");
	    roiManager("Add");
		for (i = 0; i < C; i++) {
		    selectWindow(hsTitle);
		    Stack.setPosition(i+1, izs+1, 1);
			roiManager("Select", 0);
			getStatistics(area, meanEC, minImg, max, std, histogram);
		    selectWindow(meanTrendTitle);
			setPixel(izs, i, meanEC);
		}
	}
}

// takes input 2D image inputImg and generates a 2D image titled linTrend
// with same size where the value of each pixel (x,y) 
// is equal the linear interpolation of its row
// after least square fitting I(i,j) = a(j)*i+b(j)
function computeLinearTrendAlongX(inputImg,linTrend){
	selectWindow(inputImg);
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	// fill an image with x coordinate value
	newImage("x", "32-bit", sizeX, sizeY, 1);
	for (i = 0; i < sizeX; i++) {
		for (j = 0; j < sizeY; j++) {
			setPixel(i, j, i);
		}
	}
	
	newImage("sx", "32-bit", sizeX, sizeY, 1);
	curX = 0;
	for (i = 0; i < sizeX; i++) {
		curX = curX+i;
		for (j = 0; j < sizeY; j++) {
			setPixel(i, j, curX);
		}
	}

	newImage("sx2", "32-bit", sizeX, sizeY, 1);
	curX2 = 0;
	for (i = 0; i < sizeX; i++) {
		curX2 = curX2+i*i;
		for (j = 0; j < sizeY; j++) {
			setPixel(i, j, curX2);
		}
	}
	
	newImage("sy", "32-bit", sizeX, sizeY, 1);
	curY = newArray(sizeY);
	for (i = 0; i < sizeX; i++) {
		for (j = 0; j < sizeY; j++) {
			selectWindow(inputImg);
			curY[j] = curY[j]+getPixel(i,j);
			
			selectWindow("sy");
			setPixel(i, j, curY[j]);
		}
	}
	
	newImage("sxy", "32-bit", sizeX, sizeY, 1);
	curXY = newArray(sizeY);
	for (i = 0; i < sizeX; i++) {
		for (j = 0; j < sizeY; j++) {
			selectWindow(inputImg);
			curXY[j] = curXY[j]+getPixel(i,j)*i;
			
			selectWindow("sxy");
			setPixel(i, j, curXY[j]);
		}
	}

	a = newArray(sizeY);
	b = newArray(sizeY);
	for (j = 0; j < sizeY; j++) {
		selectWindow("sxy");
		sxy = getPixel(sizeX-1,j);

		selectWindow("sx");
		sx = getPixel(sizeX-1,j);

		selectWindow("sx2");
		sx2 = getPixel(sizeX-1,j);

		selectWindow("sy");
		sy = getPixel(sizeX-1,j);

		a[j] = ( sizeX * sxy - sx*sy)/(sizeX *sx2 - sx*sx);
		b[j] = (sy - a[j] * sx)/ sizeX;
	}

	newImage(linTrend, "32-bit", sizeX, sizeY, 1);
	for (i = 0; i < sizeX; i++) {
		for (j = 0; j < sizeY; j++) {
			setPixel(i,j, i*a[j]+b[j]);	
		}
	}
	selectWindow("x");
	close();
	selectWindow("sx");
	close();
	selectWindow("sx2");
	close();
	selectWindow("sy");
	close();
	selectWindow("sxy");
	close();
}



//
function selectCorrectlySegmentedNuclei(inputImageTitle,outputImageTitle,
	avgNucleiDiameterInUm,SurfaceToVolumeRatioThresh,maxBreadth,MaxOutputObj){
	
	volRatio = 36;
	minVol = avgNucleiDiameterInUm/volRatio;
	
	selectWindow(inputImageTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	
	// collect table of object geometric features
	run("Analyze Regions 3D", "volume surface_area mean_breadth centroid"
	+" surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");

	// rename the morpholibj output table so ImageJ recognizes it as a Results table
	// the extension of the filename is not transfered to the morpholibJ table name
	if(indexOf(inputImageTitle,'.')>0){
		imgNameWOExt = inputImageTitle.substring(0,indexOf(inputImageTitle,'.'));
	}else{
		imgNameWOExt = inputImageTitle;
	}
	print(imgNameWOExt);
	
	print("sorting through "+nResults+" objects...");
	// computing distances to stack center
	print("computing distances to stack center");
	Table.rename(imgNameWOExt+"-morpho", "Results");
	
	// compute distance of each opbject centroid to the stack center
	getVoxelSize(dX, dY, dZ, u);
	iX = dX*sizeX/2;
	iY = dY*sizeY/2;
	iZ = dZ*sizeZ/2;
	
	print("iX: "+iX+" iY: "+iY+" iZ: "+iZ);
	distArray = newArray(nResults);
	for (i=0; i<nResults; i++){
	 	 oX = getResult("Centroid.X",i);
	 	 oY = getResult("Centroid.Y",i);
	 	 oZ = getResult("Centroid.Z",i);
	 	 print("oX: "+oX+" oY: "+oY+" oZ: "+oZ);
	 	 
	 	 distArray[i] = sqrt( (oX-iX)*(oX-iX) + (oY-iY)*(oY-iY) + (oZ-iZ)*(oZ-iZ) );
	 	 //curLabel = getResult("Label",i);
	 	 //print("i = "+i+" label: "+curLabel+" "+distArray[i]);
	}
	
	// sort nuclei by increasing distance to image center.
	Table.setColumn("dist", distArray);
	Table.sort("dist");
	
	// filter out nuclei that do not pass quantitative metrics
	volArray = newArray(nResults);
	surfVolRatioArray = newArray(nResults);
	breadthArray = newArray(nResults);
	labelArray = newArray(nResults);
	
	counter=0;
	labelString="";
	for (i=0; i<nResults; i++){
 	 volArray[i] = getResult("Volume",i);
 	 surfVolRatioArray[i] = getResult("SurfaceArea",i)/volArray[i];
 	 labelArray[i] = getResultString("Label",i);
 	 breadthArray[i] = getResult("MeanBreadth",i);
 	 if ((volArray[i] > minVol) 
 	 		&& (surfVolRatioArray[i]<SurfaceToVolumeRatioThresh)
 	 		&& (breadthArray[i] < maxBreadth)  
 	 		&& (counter<MaxOutputObj) ){
 	 	counter=counter+1;
		print("****");
 	 	print("Selected object # " + counter);
 	 	print("Object label "+ labelArray[i]);
 	 	print("surface to Volume ratio = "+ surfVolRatioArray[i]);
 	 	print("MeanBreadth = "+ breadthArray[i]);
 	 	print("Volue = "+ volArray[i]);
 	 	labelString=labelString+labelArray[i]+",";
 	 	print(" ");
 	 }else{
 	 	print("****");
 	 	print("Rejected object # " + counter);
 	 	print("Object label "+ labelArray[i]);
 	 	print("surface to Volume ratio = "+ surfVolRatioArray[i]);
 	 	print("MeanBreadth = "+ breadthArray[i]);
 	 	print("Volue = "+ volArray[i]);
 	 }
	}
	labelString = labelString.substring(0,labelString.length-1);
	print("Selected "+ counter +" objects with following labels: "+labelString);
	
	// generate a new image with only the selected objects
	run("Select Label(s)", "label(s)="+labelString); 
	labelImg = getTitle();
	
	run("Connected Components Labeling", "connectivity=26 type=[16 bits]");
	run("glasbey on dark");
	rename(outputImageTitle);
	selectWindow(labelImg);
	close();
	return counter;
}

// function that starts from a hyperstack containing a nuclear stain 
//(enter the index of the nucleus stain in nuclearChannel),
// 0) downsample the stack by a factor of 2 in all dimensions to speed things up
// 1) perfoms a bandpass filter to smooth nuclei
// 2) thresholds nuclei from the background 
// 3) runs an opening transformation to remove smaller objects
// 4) performs a watershed to split conjoined nuclei
// 5) performs a second opening to remove small muclei
// 6) performs a dilation to enlarge the results
// 7) restores results to original size
function segmentNuclei3D(originalImgTitle,outputImgTitle,nucleiChannel,avgNucleiDiameterInUm){
	
	// ad hoc factor for the band pass filter (threshold lengths will be :
	// (low pass) any length above <diam of a nucleus>/fftSmallFactor1
	// (hi pass) any length under <diam of a nucleus>*fftLargeFactor1
	fftSmallFactor1 = 8;
	fftLargeFactor1 = 1;
	
	// ratio of the nuclear diameter to the size of the structuring element for the first opening transformation
	opening1Factor = 10; 
	
	// objects smaller than avgNucleiVolume/minVolFactor are excluded after the first filtering step,
	// where avgNucleiVolume is computed assuming a sphere of diameter avgNucleiDiameterInUm
	minVolFactor = 4.5;
	
	// ad hoc factor for the second low pass filter (applied to the distance map
	// to be ~half the diameter of a nucleus - any length above <diam of a nucleus>/fftSmallFactor2 is retained
	fftSmallFactor2 = 2;
	
	// ratio of the nuclear diameter to the size of the structuring element for the first opening transformation
	opening2Factor = 5;
	
	// ratio of the nuclear diameter to the size of the structuring element for the last dilation
	dilationFactor = 15;
	
	selectWindow(originalImgTitle);
	run("Duplicate...", "title=hoechst duplicate channels="+nucleiChannel);
	
	//downsize image to accelerate segmentation
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	x2 = round(sizeX/2);
	y2 = round(sizeY/2);
	z2 = round(sizeZ/2);
	run("Size...", "width="+x2+" height="+y2+" depth="+z2
		+" constrain average interpolation=Bicubic");
	
	// collect um <-> voxel conversion factors
	getVoxelSize(dX, dY, dZ, u);
	
	//set limits of the bandpass smoothing filter used prior to thresholding
	fftLarge1 = round(fftLargeFactor1*avgNucleiDiameterInUm/dX);
	fftSmall1 = round((avgNucleiDiameterInUm/dX)/fftSmallFactor1);
	//print("fftSmall1: "+fftSmall1 + "; fftLarge1: "+fftLarge1);
	
	// set size of structuring element for opening to the avg nucl diameter divided
	// by opening1Factor
	openingXY1 = round(avgNucleiDiameterInUm/dX/opening1Factor);
	openingZ1 = round(avgNucleiDiameterInUm/dZ/opening1Factor);
	//print("openingXY1: "+ openingXY1+"; openingZ1 = "+ openingZ1);
	
	// set size of smoothing filter to apply to distance map before watershed
	fftLarge2 = maxOf(sizeX,sizeY)*10; // factor 10 to get arbitrary large threshold because we just want a low pass
	fftSmall2 = round(avgNucleiDiameterInUm/dX/fftSmallFactor2); // this threshold is in pixel units
	//print("fftSmall2: "+fftSmall2 + "; fftLarge2: "+fftLarge2);
	
	// set a threshold for the minimum volume for accepted objects to be the expected to 
	// nucleus volume (based on the avgNucleiDiameter) divided by the ad hoc minVolFactor
	minVolInVoxels = round(4*3.14/3* pow(avgNucleiDiameterInUm/2,3)/(dX*dY*dZ)/minVolFactor);
	//print("minimum volume in voxels: "+minVolInVoxels);
	
	// set the size of structuring element for opening to the avg nucl diameter divided
	// by opening2Factor
	openingXY2 = round(avgNucleiDiameterInUm/dX/opening2Factor);
	openingZ2 = round(avgNucleiDiameterInUm/dZ/opening2Factor);
	//print("openingXY2: "+ openingXY1+"; openingZ2 = "+ openingZ1);
	
	// set the size of the structuring element for the dilation at the end
	dilationSize = round(avgNucleiDiameterInUm/dZ/dilationFactor);
	
	//nuclei segmentation 
	
	// smooth then threshold DAPI/Hoechst channel, run an opening on the result to remove small objects
	// and aberrant links between neighnbors
	run("Bandpass Filter...", "filter_large="+fftLarge1+" filter_small="+fftSmall1+" suppress=None tolerance=5 process");
	run("Convert to Mask", "method=Default background=Dark calculate black");
	run("Fill Holes", "stack");
	run("Morphological Filters (3D)", "operation=Opening element=Ball"
		+" x-radius="+openingXY1+" y-radius="+ openingXY1+ " z-radius=" +openingZ1);
	binaryImgTitle = "filteredMasks";
	rename(binaryImgTitle);
		
	// split conjoined objects - this is done by 
	// 1) a 3D champfer distance map, followd by smoothing by FFT low-pass filter 
	// 2) an inversion of the distance map followed by rescaling to 0-255 range 
	// this makes the nuclei centers the lowest values, and background highest.
	// 3) a watershed segmentation of the inverted map (limited to segmented nuclei)
	run("Chamfer Distance Map 3D", 
		"distances=[Quasi-Euclidean (1,1.41,1.73)] output=[16 bits] normalize");
	run("Bandpass Filter...", "filter_large="+fftLarge2+" filter_small="+fftSmall2+
		" suppress=None tolerance=5 process");
	run("Invert", "stack");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setMinAndMax(min, max);
	run("8-bit");
	invertedImgTitle = getTitle();
	run("Classic Watershed", "input="+invertedImgTitle+" mask="+binaryImgTitle+" use min=0 max=255");
	
	// filter out objects too small to be nuclei
	selectWindow("watershed");
	run("Label Size Filtering", "operation=Greater_Than size="+minVolInVoxels);
	run("glasbey on dark");
	
	// run an opening filter in order to smooth the shape of the nuclei
	selectWindow("watershed-sizeFilt");
	run("Morphological Filters (3D)", "operation=Opening element=Ball"
	+" x-radius="+openingXY2+" y-radius="+openingXY2+" z-radius="+openingZ2);
	
	// dilate the nuclei to their full size
	run("Morphological Filters (3D)", "operation=Dilation element=Ball"
	+" x-radius="+3+" y-radius="+3+" z-radius="+0);
	
	// rescale the segmentation back to its original size
	run("Size...", "width="+sizeX+" height="+sizeY+" depth="+sizeZ+" interpolation=None");
	setMinAndMax(0, 65535);
	run("16-bit");
	res = getTitle();
	
	// close intermediates
	selectWindow("hoechst");
	close();
	selectWindow("filteredMasks");
	close();
	selectWindow("filteredMasks-dist");
	close();
	selectWindow("watershed");
	close();
	selectWindow("watershed-sizeFilt");
	close();
	selectWindow("watershed-sizeFilt-Opening");
	close();
	
	// output segmentation result
	Stack.setPosition(2, 1, 1);
	setMinAndMax(0, 20);
	rename(outputImgTitle);
}

