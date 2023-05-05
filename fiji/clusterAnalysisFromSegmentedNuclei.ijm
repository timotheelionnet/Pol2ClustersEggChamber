c// this macro detects clusters from a dataset of ovary egg chamber z-stacks where the nuclei are already segmented.
// these nuclei segmentations can be generated from the raw z-stacks using the macro segmentNuclei.ijm

// the expected input to this macro is a data folder which should have the following architecture:
// condition names and sample names are arbitrary, but the name of the sample folder should match
// the begining of the name of the file holding the z-stack in the effChanberTIF folder.

// <dataFolder/> 
//		|__ <condition1/>
//				|__<sample1/>
//						|__ <eggChamberTIF/>	
//								|__sample1ZcorrFinalNucMask.tif			
//				|__<sample2/>
//						|__ <eggChamberTIF/>	
//								|__sample2ZcorrFinalNucMask.tif	
//				|__ ...
//
//		|__ <condition2/>
//				|__<sample1/>
//						|__ <eggChamberTIF/>	
//								|__sample1ZcorrFinalNucMask.tif		
// 		|__ ...

// Each z-stack is expected to be a multi-channel hyperstack, with n channels of data, and the *last* channel 
// should hold the result of the nuclei segmentations, ie. background set to zero and
// each segmented nucleus set to its unique integer ID value.  

// The macro goes through the list of conditions and samples, and
// - loads each egg chamber z-stack (the image "sample1ZcorrFinalNucMask.tif")
// - generates a cropped z-stack around each individual nucleus.
// - within each segmented nucleus, segments the nucleoplasm based on the input channel set as pol2 
//		(keeping only the nucleoplasm helps eliminate nucleoli which tend to have much lower level of Pol2 
//		and thus make cluster detection more difficult).
// - detects clusters within the nucleoplasm mask based on the input channel set as pol2.
// - measures statistics of the clusters in all channels (size, intensity, etc)
// - saves the statistics as csv files for each channel and each nucleus:
//			- e.g. for nucleus 1, in channel 1:
//				- C1_nuc1_plasmInt.csv stores geometry and intensity data about the nucleoplasm (single row table)
//				- C1_nuc1_clustInt_raw.csv stores geometry and raw intensity data 
//						about the clusters found in this nucleus
//				- C1_nuc1_clustInt_plasmCorr.csv stores geometry and intensity data 
//						(from which the average nucleoplasm intensity level has been subtracted out) 
//						 about the clusters found in this nucleus
// - saves each cropped nucleus stack into a single z-stack, along with the segmentation results
//		e.g. for nucleus 1: nuc1_masks.tif in the folder nucTIF. 
//		extra channels are appended at the end of the hyperstack (in this order) for:
//			- the nucleus segmentation (copied from the input); 
//			- the nucleoplasm segmentation
//			- the cluster segmentation
// - also saves the same z-stack, from which the average nucleoplasm intensity has been subtracted in each channel
// 		under the name nuc1_masks_plasmCorr.tif (for nucleus 1).
//		the plasmCorr z-stack is a 32 bit format since many pixels end up negative upon
//		background subtraction.


// The end result is something like:
// <dataFolder/> 
//		|__ <condition1/>
//				|__<sample1/>
//						|__ <eggChamberTIF/>	
//								|__sample1ZcorrFinalNucMask.tif
//						|__ <nucTIF/>
//								|__ nuc1_masks.tif	
//								|__ nuc1_masks_plasmCorr.tif
//								|__ nuc2_masks.tif	
//								|__ nuc2_masks_plasmCorr.tif
//								|__ ...
//						|__ <nucCSV/>	
//								|__ C1_nuc1_clustInt_plasmCorr.csv
//								|__ C1_nuc1_clustInt_raw.csv
//								|__ C1_nuc1_plasmInt.csv
//								|__ C2_nuc1_clustInt_plasmCorr.csv
//								|__ C2_nuc1_clustInt_raw.csv
//								|__ C2_nuc1_plasmInt.csv
//								|__ C3_nuc1_clustInt_plasmCorr.csv
//								|__ C3_nuc1_clustInt_raw.csv
//								|__ C3_nuc1_plasmInt.csv
//								|__ ...
//								|__ C1_nuc2_clustInt_plasmCorr.csv
//								|__ C2_nuc2_clustInt_plasmCorr.csv
//								|__ C3_nuc2_clustInt_plasmCorr.csv
//								|__ ...
						
//				|__<sample2/>
//						|__ ...

//				|__ ...
//
//		|__ <condition2/>
//				|__<sample1/>
//						|__ ...
//				|__ ...
	
// 		|__ ...

var border = 1; // padding (in pixels) around each nuclei in its individually cropped stack
var saveStacks = 0; // harcoded variable; suggested setting 0; set to 1 to save stacks with only the nucleus segmentation added (these stacks are generated by cropNuclei as an intermediate step redundant with the final output)
var pol2Channel = 4;
var MedianFilterRadiusClusters = 2;
var MedianFilterRadiusBackground = 50;
var threshFactor = 3;
var verbose = 1;

macro "clusterAnalysisFromSegmentedNuclei" {
	setBatchMode(true);
	inputParameters();
	
	segmentationSubFolderName = "eggChamberTIF/";
	segmentationFileSuffix = "ZcorrFinalNucMask.tif";
	csvSaveDirName = "nucCSV/"; 
	tifSaveDirName = "nucTIF/";
	
	// define input folder (should be the output folder of the nucleus segmentation step
	// performed by segmentNuclei). Results will be saved in the same root folder
	//inFolder = "/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2/";
	inFolder = getDirectory("choose the input directory");
	
	// loop through conditions (first level of subfolders from input folder)
	conditionList = getSubdirList(inFolder);
	print("Found the following conditions:");
	for(i=0;i<conditionList.length;i++){
		print(conditionList[i]);
	}
	
	for(i=0;i<conditionList.length;i++){
		// loop though eggchambers (second level of subfolders from input folder)
		sampleList = getSubdirList(inFolder+conditionList[i]);
		
		for(j=0;j<sampleList.length;j++){
			
			// reconstruct file dir
			curImgFilePath = inFolder + conditionList[i] + sampleList[j] + segmentationSubFolderName;
			// reconstruct file name
			// remove last character which is a / because sampleList[j] is a directory
			s = substring(sampleList[j], 0,lengthOf(sampleList[j])-1); 
			curImgFileName = s + "ZcorrFinalNucMask.tif";	
			
			if( File.exists(curImgFilePath+curImgFileName)){
				// load eggchamber nucleus segmentation file
				print("cond "+i+"/sample "+j+"; loading "+curImgFilePath+curImgFileName +" ...");
				open(curImgFilePath+curImgFileName);
				rename("mergedSegmentation");
				
				// generate a cropped image for each nucleus that holds the original image channels 
				// plus the nuclear mask as an extra color channel
				// images are called "croppedNucleus"+label+"wSeg_masks.tif"
				print("cropping nuclei...");
				nucleiFound = cropNuclei("mergedSegmentation",border,saveStacks);
				
				// loop through nuclei, segment nucleoplasm and clusters
				for (k = 0; k < nucleiFound.length; k++) {
					label = nucleiFound[k];
					print("processing nucleus "+label+"...");
					curImg = "croppedNucleus"+label+"wSeg";
					getDimensions(originalImgSizeX, originalImgSizeY, segmentationC, originalImgSizeZ, originalImgF);
					maskChannel = segmentationC;
					
					// compute nucleus geometry in cropped image - the centroid coordinates can 
					// then be used to convert local coordiantes (i.e. relative to the cropped img)
					// into global ones (i.e. relative to the global egg chamber img).
					csvSaveDir = inFolder + conditionList[i] + sampleList[j] + csvSaveDirName;
					if(!File.exists(csvSaveDir)){
						File.makeDirectory(csvSaveDir);
					}
					nucleiGeometryMeasurement(curImg,maskChannel,csvSaveDir+"nuc"+label+"_localGeom.csv"); 
					
					// segment nucleoplasm - output is an 8-bit mask image called maskedNucleoplasmImgTitle
					// where background is 0 and nucleoplasm is 255.
					print("segmenting nucleoplasm of nucleus "+label+"...");
					maskedNucleoplasmImgTitle = "maskedPlasm";
				  	segmentNucleoplasm(curImg,pol2Channel,maskedNucleoplasmImgTitle);
			
				  	//subtract the average signal value in the nucleoplasm from all channels
				  	print("subtracting background of nucleus "+label+"...");
				  	imgBgCorr = "backgroundCorrectedImg";
					subtractBackground32bit3D(curImg,maskedNucleoplasmImgTitle,imgBgCorr);
					
				  	// segment clusters - output is a 16 bit image called maskedClustersImgTitle
				  	// where background is 0 and each cluster carries an individual ID
				  	print("detecting clusters of nucleus "+label+"...");
				  	maskedClustersImgTitle = "clusterMask";
				  	//detectClusters(imgBgCorr,maskedClustersImgTitle,pol2Channel,
				  	//	MedianFilterRadiusClusters,MedianFilterRadiusBackground,
				  	//	maskedNucleoplasmImgTitle,threshFactor);
				  	detectClusters3D(imgBgCorr,maskedClustersImgTitle,maskChannel,pol2Channel,
						MedianFilterRadiusClusters,MedianFilterRadiusBackground,
						maskedNucleoplasmImgTitle,threshFactor);
					
					// compute parameters of clusters from raw img and save results
				  	print("computing cluster stats of nucleus "+label+" from raw stack...");
				  	computeClusterStats(curImg,maskedClustersImgTitle,csvSaveDir,label,"_clustInt_raw.csv");
					
				  	// compute parameters of clusters from nucleoplasm subtracted image and save results
				  	print("computing cluster stats of nucleus "+label+" from nucleoplasm-corrected stack...");
				  	computeClusterStats(imgBgCorr,maskedClustersImgTitle,csvSaveDir,label,"_clustInt_plasmCorr.csv");
			
				  	// generate output image that includes the original hyperstack plus 3 extra channels appended at the end:
				  	// (1) the masks of the nucleus,
				  	// (2) the mask of the nucleoplasm and 
				  	// (3) the mask of the Pol II clusters	
				  	print("generating output image of nucleus "+label+"...");
				  	imgOutName = "nuc"+label +"_masks.tif";
				  	selectWindow(maskedNucleoplasmImgTitle);
				  	run("16-bit");
				  	addChannelToImg(curImg,maskedNucleoplasmImgTitle,"tmp",1);
				  	addChannelToImg("tmp",maskedClustersImgTitle,imgOutName,1);
					selectWindow("tmp");
					close();
					selectWindow(curImg);
					close();
					
				  	// collect intensity measurements on the nucleoplasm
					nucleoplasmIntensityAndGeometryMeasurementAllChannels(
						csvSaveDir,"nuc"+label,imgOutName,maskChannel+1);
					
					// save image with masks
					selectWindow(imgOutName);
					tifSaveDir = inFolder + conditionList[i] + sampleList[j] + tifSaveDirName;
					if(!File.exists(tifSaveDir)){
						File.makeDirectory(tifSaveDir);
					}
				  	save(tifSaveDir+imgOutName);
				  	close(imgOutName);
			
					// save nucleoplasm-background corrected image
					bgCorrImgOutName = "nuc"+label +"_masks_plasmCorr.tif";
				  	print("generating output background-corrected image of nucleus "+label+"...");
				  	addChannelToImg(imgBgCorr,maskedNucleoplasmImgTitle,"tmp",0);
				  	addChannelToImg("tmp",maskedClustersImgTitle,bgCorrImgOutName,0);
					selectWindow(bgCorrImgOutName);
				  	save(tifSaveDir+bgCorrImgOutName);
				  	close(bgCorrImgOutName);	
				  	
				  	print("saving parameters to file...");
					saveClusterParametersToLogFile(curImgFilePath+curImgFileName,
						inFolder+ conditionList[i] + sampleList[j] + "clusterAnalysisParameters.txt",
						originalImgSizeX, originalImgSizeY, segmentationC - 1, originalImgSizeZ,
						border,nucleiFound,
						pol2Channel,threshFactor,MedianFilterRadiusClusters,
						MedianFilterRadiusBackground);
				}
				close("mergedSegmentation");
			}else{
				print("Could not find file "+curImgFilePath+curImgFileName+"; skipping...");
			}
		}
	}			
	setBatchMode("exit and display");
  	print("done.");			
}


// self explanatory
function inputParameters(){
	Dialog.create("Cluster Detection Parameters");
	
	Dialog.addNumber("Pol2 color channel:", pol2Channel);	
	Dialog.addNumber("Border around nuclei in cropped stacks (pixels): ", border);	
	Dialog.addNumber("Radius of Median Filter to smooth cluster signal: ", MedianFilterRadiusClusters);
	Dialog.addNumber("Radius of Median Filter to average out background: ", MedianFilterRadiusBackground);
	Dialog.addNumber("Threshold Factor to call clusters (multiples of nucleoplasm StdDev): ",threshFactor);
	Dialog.show();
	
	pol2Channel = Dialog.getNumber();
	border = Dialog.getNumber();
	MedianFilterRadiusClusters = Dialog.getNumber();
	MedianFilterRadiusBackground =  Dialog.getNumber();
	threshFactor =  Dialog.getNumber();

	print("*******************************");
	print("Pol II color channel: "+pol2Channel);
	print("Border around nuclei in cropped stacks (pixels): ",border);
	print("Radius of Median Filter to smooth cluster signal: ",MedianFilterRadiusClusters);
	print("Radius of Median Filter to average out background: ",MedianFilterRadiusBackground);
	print("Threshold Factor to call clusters (multiples of nucleoplasm StdDev): ",threshFactor);
}
function detectClusters3D(curImg,maskedClustersImgTitle,nucMaskChannel,pol2channel,
	MedFiltRadiusCluster,MedFiltRadiusBg,maskedNucleoplasmImgTitle,threshF){
	
	// median filter Pol II channel for cluster detection
	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+pol2Channel);
	run("Median...", "radius="+MedFiltRadiusCluster+" stack");
	rename("pol2median");

	run("Intensity Measurements 2D/3D", 
		"input=pol2median labels="+maskedNucleoplasmImgTitle+" mean stddev");

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
	
	//isolate clusters
	selectWindow("pol2median");
	run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
	run("glasbey on dark");
	rename("preSegmentedClusters");
	close("pol2median");
	
	//copy nucleus mask channel
	selectWindow(curImg);
	run("Duplicate...", "duplicate channels="+nucMaskChannel);
	rename("tmpNucMasks");
	
	// measure max intensity in each cluster in the hoechst channel
	run("Intensity Measurements 2D/3D", "input=tmpNucMasks labels=preSegmentedClusters max");
	close("tmpNucMasks");
	
	// rename the morpholibj output table so ImageJ recognizes it as a Results table
	// the extension of the filename is not transfered to the morpholibJ table name
	Table.rename("tmpNucMasks-intensity-measurements", "Results");
	
	// find all clusters with non-zero max intensity in the nucleus mask channel and collect their
	// indices in a string to isolate them later
	labelString = "label(s)=";
	ctr = 0;
	for (i=0; i<nResults; i++){
		if(getResult("Max",i)>0){
			labelToKeep = getResultString("Label",i);
			if(ctr == 0){
				labelString = labelString+labelToKeep;
			}else {
				labelString = labelString+","+labelToKeep;
			}
			ctr = ctr+1;
		}
	}
	close("Results");
	
	// keep only clusters with non-zero max intensity 
	//(i.e. clusters that overlap w the nucleus mask at least partially)
	selectWindow("preSegmentedClusters");
	if (labelString != "label(s)="){
		run("Select Label(s)", labelString);
		print("Selected "+ ctr +" clusters with following labels: "+labelString);
	}else{
		print("Found no clusters to include.");
	}
	rename(maskedClustersImgTitle);
	close("preSegmentedClusters");
}

function saveClusterParametersToLogFile(originalImgTitle,savePath,
	sizeX, sizeY, C, sizeZ,
	border,nucleiFound,
	pol2Channel,threshFactor,MedianFilterRadiusClusters,
	MedianFilterRadiusBackground){

	f = File.open(savePath);
   	// use d2s() function (double to string) to specify decimal places 
   	//for (i=0; i<=2*PI; i+=0.01)
    //  print(f, d2s(i,6) + "  \t" + d2s(sin(i),6) + " \t" + d2s(cos(i),6));
    print(f, "****** Input image: \r\n");
    print(f,"file: "+ originalImgTitle + "\r\n");
    print(f,"dimensions : "+ sizeX +" (X), "+sizeY+" (Y), "+ sizeZ + " (Z), " +C+ " (Channels).\r\n");
    print(f," \r\n");
    
    print(f,"****** Nucleus Detection Settings:\r\n");
	print(f, "border: "+ border+ "\r\n");
	print(f, "nucleiFound: "+ nucleiFound.length+ "\r\n");
	print(f," \r\n");
	
	print(f,"****** Cluster Detection Settings:\r\n");
	print(f, "pol2Channel: "+ pol2Channel+ "\r\n");
	print(f, "threshFactor: "+ threshFactor+ "\r\n");
	print(f, "MedianFilterRadiusClusters: "+ MedianFilterRadiusClusters+ "\r\n");
	print(f, "MedianFilterRadiusBackground: "+ MedianFilterRadiusBackground+ "\r\n");
	print(f," \r\n");
	File.close(f);
	
}


// takes an image (inputWindowName), duplicates the channel with nucleoplasm mask, 
// then measures 1) intensity etc and 2) geometry on all Channels
// and saves the resulting tables in the folder savePath using names like 
// C1_rootName_plasmInt.csv, C2_rootName_plasmInt.csv, etc for intensity features
// single file rootName-nucleoplasmGeom.csv for geometry features 
function nucleoplasmIntensityAndGeometryMeasurementAllChannels(savePath,rootName,inputWindowName,maskChannel){
	selectWindow(inputWindowName);
	run("Duplicate...", "duplicate channels="+maskChannel);
	rename("tmpDuplicate1");

	run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity"+
			" euler_number bounding_box centroid equivalent_ellipsoid ellipsoid_elongations"+
			" max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");

	saveAs("Results",savePath+rootName+"_plasmGeom.csv");
	run("Close");
	
	selectWindow(inputWindowName);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	for (i = 1; i <= C; i++) {
	    selectWindow(inputWindowName);
	    run("Duplicate...", "duplicate channels="+i);
		rename("tmpDuplicate2");
		run("Intensity Measurements 2D/3D", "input=tmpDuplicate2 labels=tmpDuplicate1"+ 
			" mean stddev max min median mode skewness kurtosis numberofvoxels volume"+
			" neighborsmean neighborsstddev neighborsmax neighborsmin neighborsmedian"+
			" neighborsmode neighborsskewness neighborskurtosis");
		
		saveAs("tmpDuplicate2-intensity-measurements",savePath+"C"+i+"_"+rootName+"_plasmInt.csv");
		run("Close");
		selectWindow("tmpDuplicate2");
		close();
	}
	
	print("Nucleoplasm Intensity and Geometry results saved");
	selectWindow("tmpDuplicate1");
	close();
}

function computeClusterStats(dataImg,labelImg,saveDir,nucleusNumber,csvSuffix){
	
	selectWindow(dataImg);
	getDimensions(w, h, c, nzs, f);
	for (i = 1; i <= c; i++) {
		selectWindow(dataImg);
		run("Duplicate...", "duplicate channels="+i);
		rename("tmpChannel");
		run("Intensity Measurements 2D/3D", "input=tmpChannel labels="
			+labelImg+" mean stddev max min median mode skewness kurtosis numberofvoxels volume");
		saveAs("Results", saveDir+"C"+i+"_nuc"+nucleusNumber +csvSuffix);
		run("Close");
		selectWindow("tmpChannel");
		close();
	}
}

// from a hypertack curImg, containing a channel with Pol II data (pol2channel),
// detects Pol II clusters as follows and saves a mask stack of the detected clusters under the title maskedClustersImgTitle.
// cluster detection proceeds as follows:
// 1) filter out salt and pepper noise in Pol II channel (median filtering with radius MedFiltRadiusCluster)
// 2) subtract the nuclear background from the filtered Pol II channel (this requires a stack holding a mask of the nucleoplasm as an input, maskedNucleoplasmImgTitle)
// 3) detect clusters as regions where intensity is threshF times the std of the intensity of the background corrected image.
// #2) and #3) are performed slice by slice to accound for loss of signal with increasing depth.
// Because nurse cells nuclei have the background is computed by:
// 2.1 creating image that is 0 outside nucleoplasm, smoothed pol2 signal inside
// 2.2 fill in the holes of the image in 2.1 with the average intensity in the nucleoplasm
function detectClusters(curImg,maskedClustersImgTitle,pol2channel,MedFiltRadiusCluster,MedFiltRadiusBg,maskedNucleoplasmImgTitle,threshF){

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
			close("backgroundSubstracted");
			close("tempPol2");
			close("filledNucleoplasm");
			close("medFilteredBackground");
			close("outsideNucleoplasm");
		}
	}
	selectWindow("clusterMasks");
	// run connected components labeling so that each nucleus gets an individual ID
	run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
	run("glasbey on dark");
	rename(maskedClustersImgTitle);
	
	// close intermediate
	close("clusterMasks");
	close("plasmMask");
	close("bgCorrImg");
	close("pol2median");
}

// takes curImg hyperstack as input, using the nucleoplasm mask stack
// at each slice it subtractes the average nucleoplasm value to the whole plane
function subtractBackground32bit3D(curImg,maskedNucleoplasmImgTitle,imgBgCorr){
	selectWindow(curImg);
	getDimensions(w, h, c, nzs, f);
	
	// compute mean in the nucleoplasm mask in each channel
	plasmMean = newArray(c);
	for (ic = 1; ic <= c; ic++) {
		run("Duplicate...", "duplicate channels="+ic);
		rename("curChannel");
		run("Intensity Measurements 2D/3D", "input=curChannel labels="+maskedNucleoplasmImgTitle+" mean");

		// rename the morpholibj output table so ImageJ recognizes it as a Results table
		// the extension of the filename is not transfered to the morpholibJ table name
		Table.rename("curChannel-intensity-measurements", "Results");
	
		// extract  mean
		plasmMean[ic-1] = getResult("Mean",0);
		close("curChannel");
		close("Results");
	}
	
	//duplicate the image
	selectWindow(curImg);
	run("Duplicate...", "duplicate");
	rename(imgBgCorr);
	run("32-bit");
	
	// subtract background from duplicate in each channel
	// only goes to c-1 because the last channel is the nucleus segmentation result
	selectWindow(imgBgCorr);
	run("Select All");
	for(ic=0;ic<c-1;ic++){
		for (izs = 0; izs < nzs; izs++) {
			Stack.setPosition(ic+1, izs+1, 1);
			run("Subtract...", "value="+plasmMean[ic]+" slice");
		}
	}
}

// starts from an image curImg where the Pol II signal is in channel pol2ChannelNum
// threshold the nuclei based on Pol II signal, then
// duplicate last channel holding nuclei IDs
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

function nucleiGeometryMeasurement(inputWindowName,maskChannel,saveFileName){
	selectWindow(inputWindowName);
	run("Duplicate...", "duplicate channels="+maskChannel);
	rename("tmpDuplicate1");

	run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity"+
			" euler_number bounding_box centroid equivalent_ellipsoid ellipsoid_elongations"+
			" max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");

	saveAs("Results",saveFileName);
	run("Close");
	selectWindow("tmpDuplicate1");
	close();
}

// this function loops through all nuclei objects (defined by masks in the last color channel of the hyperstack inputImgTitle)
// and generates one stack per nucleus (called "croppedNucleus"+i+"wSeg" for nucleus i)
// that retains all the color channels of the input stack; 
// the last channel only contains the mask corresponding to the nucleus of interest (other nuclei masks are erased)
// border is the padding in pixels around each cropped nucleus
function cropNuclei(inputImgTitle,border,saveStacks){
	selectWindow(inputImgTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	segmentationChannel = C;
	run("Duplicate...", "duplicate channels="+C);
	rename("duplicateNucleiID");

	selectWindow(inputImgTitle);
	C2 = C-1;
	run("Duplicate...", "duplicate channels=1-"+C2);
	rename("dataStack");
	
	selectWindow("duplicateNucleiID");
	run("Z Project...", "projection=[Max Intensity]");
	getStatistics(area, mean, min, max, std);
	maxNucleiFound = max;
	nIdx = 0;
	nucleiFound = newArray(1000);
	for (i = 1; i < maxNucleiFound+1; i++) {
		if (isNucleusIndexPresent("duplicateNucleiID",i)==1){
			nucleiFound[nIdx] = i;
			nIdx = nIdx+1;
			getLabelBoundingBox3("duplicateNucleiID","dataStack","croppedNucleus"+i+"wSeg",i,border);
	
			// save image with masks
			selectWindow("croppedNucleus"+i+"wSeg");
			if (saveStacks == 1){
				save(saveDir+"croppedNucleus"+i+"wSeg_masks.tif");	
			}
		}
	}
	nucleiFound = Array.trim(nucleiFound,nIdx);
	
	// close intermediate
	selectWindow("duplicateNucleiID");
	close();
	selectWindow("dataStack");
	close();
	selectWindow("MAX_duplicateNucleiID");
	close();
	print("number of nuclei =" + nIdx);
	return nucleiFound;
	
}

function isNucleusIndexPresent(imgName,idx){
	selectWindow(imgName);
	run("Duplicate...", "duplicate");
	rename("duplicateIdx");
	setAutoThreshold("Default dark");
	setThreshold(idx, idx, "raw");
	run("Convert to Mask", "method=Default background=Dark black");
	run("Z Project...", "projection=[Max Intensity]");
	rename("maxDuplicateIdx");
	getStatistics(area, mean, min, max, std, histogram);
	close("duplicateIdx");
	close("maxDuplicateIdx");
	if(max == 255){
		return 1;
	}else {
		return 0;
	}
}

// finds the rectangular bounding box around a 3D object labeled with the value label in the input stack inputObjectIDsTitle
// generates a rectangular ROI around the object including a padding radius of border pixels around it.
// crops the data stack inputImgDataTitle accroding to that ROI
function getLabelBoundingBox3(inputObjectIDsTitle,inputImgDataTitle,outputImgTitle,label,border){
	selectWindow(inputObjectIDsTitle);
	run("Select All");
	
	// duplicate nuclei ID image
	run("Duplicate...", "duplicate");
	
	// threshold stack to keep only current nucleus
	setThreshold(label, label);
	run("Convert to Mask", "method=Default background=Dark black");
	rename("isolatedNucleus"+label);
	
	// image size
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	// find the lowest plane holding the nucleus (zmin) and the highest (zmax)
	zmin = sizeZ;
	zmax = 0;
	for (z = 1; z <= sizeZ; z++)
	{
		Stack.setPosition(1, z, 1);
		getStatistics(area, mean, min, max, std, histogram);
		if(max == 255){
			if (zmin>z) {
				zmin=z;
			}
			if (zmax<z) {
				zmax=z;
			}
		}
	}
	print("zmin = "+zmin+"; zmax = "+zmax);

	// add border
	if((zmin-border)<1){
		zmin = 1;
	}else {
		zmin = zmin-border;
	}

	if((zmax+border)>sizeZ){
		zmax = sizeZ;
	}else {
		zmax = zmax+border;
	}
	print("(with border) zmin = "+zmin+"; zmax = "+zmax);
	
	// find xy bounding box
	run("Z Project...", "projection=[Max Intensity]");
	// add padding in xy by running a series of border dilations
	for (i = 0; i < border; i++) {
		run("Dilate");
	}
	roiManager("reset");
	setThreshold(100, 255);
	run("Create Selection");
	run("To Bounding Box");
	roiManager("Add");

	// select the original data stack and duplicate a cropped version
	selectWindow(inputImgDataTitle);
	roiManager("Select", 0);
	run("Duplicate...", "duplicate slices="+zmin+"-"+zmax);
	rename("croppedNucleus"+label);

	// cropped image size
	getDimensions(w, h, c, z, f);
	
	//convert the thresholded nucleus image back to its original value
	selectWindow("isolatedNucleus"+label);
	run("16-bit");
	run("Multiply...", "value="+label+" stack");
	run("Divide...", "value=255.000 stack");
	roiManager("Select", 0);
	run("Duplicate...", "duplicate range="+zmin+"-"+zmax);
	rename("isolatedCroppedNucleus"+label);
	
	// create cropped image holding initial data + nucleus segmentation
	addChannelToImg("croppedNucleus"+label,"isolatedCroppedNucleus"+label,outputImgTitle,0);

	// close intermediates
	selectWindow("isolatedNucleus"+label);
	close();
	selectWindow("MAX_isolatedNucleus"+label);
	close();
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

// collect names of subfolders within inFolder and
// store into subDirList, so that the path to subFolder i is
// inFolder+subDirList[i].
function getSubdirList(inFolder){
	dirList = getFileList(inFolder);
	subDirList = newArray(1000);
	ctr = 0;
	print(" ");
	for (i = 0; i < lengthOf(dirList); i++) {
	    curSubDirList = getFileList(inFolder+dirList[i]);
    	if (curSubDirList.length>0){
    		subDirList[ctr] = dirList[i];
			ctr = ctr+1;
    	}
	}
	subDirList = Array.trim(subDirList, ctr);
	return subDirList;
}

