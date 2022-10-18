macro "test Nuclei Segmentation"{
	setBatchMode(true);
	hoechstChannel = 1;
	fftSmall = 60; // minimum size of objects retained by the bandpass filter used to segment nuclei
	fftLarge = 2000; // max size
	
	nMaxObjectVolume = 15;
	SurfaceToVolumeRatioThreshold = 1.1;
	maxBreadth = 30;
	minVolume = 200;
	MaxOutputObj = 15;
	
	print("segmenting nuclei...");
	originalImgTitle = getTitle();	
	segmentNuclei3D(originalImgTitle,"firstSegResult",hoechstChannel,fftSmall,fftLarge,minVolume);
	
	print("isolating nurse cells from follicle cells");
	selectCorrectlySegmentedNuclei("firstSegResult","cleanNuclei",
		minVolume,SurfaceToVolumeRatioThreshold,maxBreadth,MaxOutputObj);	
	
	setBatchMode("exit and display");
  	print("done.");
}

function selectCorrectlySegmentedNuclei(inputImageTitle,outputImageTitle,
	minVol,SurfaceToVolumeRatioThresh,maxBreadth,MaxOutputObj){
	
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

function segmentNuclei3D(originalImgTitle,outputImgTitle,nucleiChannel,fftSmall,fftLarge,minVol){
	// duplicate hoechst channel
	selectWindow(originalImgTitle);
	run("Duplicate...", "title=hoechst duplicate channels="+nucleiChannel);
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	
	//downsize image to accelerate segmentation
	x2 = round(sizeX/2);
	y2 = round(sizeY/2);
	z2 = round(sizeZ/2);
	
	run("Size...", "width="+x2+" height="+y2+" depth="+z2+" constrain average interpolation=None");
	
	// threshold DAPI/Hoechst channel
	run("Make Binary", "method=Shanbhag background=Dark black");
	run("Fill Holes", "stack");
	//run("Morphological Filters (3D)", "operation=Erosion element=Ball x-radius=6 y-radius=6 z-radius=4");
	run("Morphological Filters (3D)", "operation=Opening element=Ball x-radius=4 y-radius=4 z-radius=2");
	binaryImgTitle = "filteredMasks";
	rename(binaryImgTitle);
		
	// compute distance map to do a watershed segmentaion of conjoined objects
	run("Chamfer Distance Map 3D", "distances=[Quasi-Euclidean (1,1.41,1.73)] output=[16 bits] normalize");
	run("Bandpass Filter...", "filter_large="+fftLarge+" filter_small="+fftSmall+" suppress=None tolerance=5 process");
	run("Invert", "stack");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setMinAndMax(min, max);
	run("8-bit");
	invertedImgTitle = getTitle();
	run("Classic Watershed", "input="+invertedImgTitle+" mask="+binaryImgTitle+" use min=0 max=255");
	
	// remove small objects
	selectWindow("watershed");
	run("Label Size Filtering", "operation=Greater_Than size="+minVol);
	rename("newWatershed");
	selectWindow("watershed");
	close();
	selectWindow("newWatershed");
	rename("watershed");
	
	//run a few dilation/opening operations to smooth the volume of the objects.
	run("Morphological Filters (3D)", "operation=Dilation element=Ball x-radius=5 y-radius=5 z-radius=4");
	rename("newWatershed");
	selectWindow("watershed");
	close();
	selectWindow("newWatershed");
	rename("watershed");
	run("Morphological Filters (3D)", "operation=Opening element=Ball x-radius=12 y-radius=12 z-radius=2");
	rename("newWatershed");
	selectWindow("watershed");
	close();
	selectWindow("newWatershed");
	rename("watershed");
	run("Morphological Filters (3D)", "operation=Dilation element=Ball x-radius=4 y-radius=4 z-radius=0");
	rename("newWatershed");
	selectWindow("watershed");
	close();
	selectWindow("newWatershed");
	rename("watershed");
	
	// remove nuclei that are cropped by the image border
	run("Remove Border Labels", "left right top bottom");
	
	// return image to size
	run("Size...", "width="+sizeX+" height="+sizeY+" depth="+sizeZ+" interpolation=None");
	run("glasbey on dark");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setMinAndMax(min, max);
	run("16-bit");
	rename(outputImgTitle);
	
	// close intermediates
	selectWindow("hoechst");
	close();
	selectWindow("watershed");
	close();
	selectWindow(binaryImgTitle);
	close();
	selectWindow(binaryImgTitle+"-dist");
	close();
}