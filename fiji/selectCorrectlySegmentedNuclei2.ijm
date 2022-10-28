// in order to sort out nurse cells nuclei (huge) from follicle cells nuclei (small)
// and to avoid grossly missegmented objects, this function loops through the largest nMaxObjectVolume nuclei 
// and selects the largest MaxOutputObject ones that have a surface:volume ratio lower than SurfaceToVolumeRatioThreshold
// nMaxObjectVolume (default 20), MaxOutputObject (15), and SurfaceToVolumeRatioThreshold (1.2) should be global variables
// that are required by the function.
function selectCorrectlySegmentedNuclei2(inputImgTitle,maskImgTitle,
	hoechstChannel,outputImgTitle
	minVolume,maxVolume, maxSurfToVolRatio,minSphericity,maxLabel,
	minMeanBreadth,maxMeanBreadth, maxSkewness,maxKurtosis,
	minCV,maxCV,useDefault){
	
	if(useDefault == 1){
		maxVolume = 5000;
		minVolume = 100; 	
		maxSurfToVolRatio = 1; 	
		minSphericity = 0.4;
		//maxLabel = 15; //17 for less stringent
		//maxLabel = 17; //this is the best, 15 for more stringent
		maxLabel = 50; //essentially not a criterion, use for higher recall, lower precision
		minMeanBreadth = 7;
		maxMeanBreadth = 25;
		maxSkewness = 1;
		maxKurtosis = 3;
		maxCV = 0.45;
		minCV = 0.15;
	}
	
	// remove masks that touch any lateral boundary
	selectWindow(maskImgTitle);
	run("Remove Border Labels", "left right top bottom");
	rename(maskImgTitle);

	// collect table of object geometric features
	run("Analyze Regions 3D", "volume surface_area mean_breadth sphericity euler_number" 
	+" bounding_box centroid equivalent_ellipsoid ellipsoid_elongations max._inscribed"
	+" surface_area_method=[Crofton (13 dirs.)] euler_connectivity=26");

	// rename the morpholibj output table so ImageJ recognizes it as a Results table
	// the extension of the filename is not transfered to the morpholibJ table name
	if(indexOf(inputImgTitle,'.')>0){
		imgNameWOExt = inputImgTitle.substring(0,indexOf(inputImgTitle,'.'));
	}else{
		imgNameWOExt = inputImgTitle;
	}
	print(imgNameWOExt);
	Table.rename(imgNameWOExt+"-morpho", "Results");
	
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
		l = getResult("Label",i);
		if(l > maxLabel){
			keepCurrent = 0;
		}
		mb = getResult("MeanBreadth",i);
		if((mb<minMeanBreadth) || (mb>maxMeanBreadth)){
			keepCurrent = 0;
		}
		
		if(keepCurrent == 1){
			labelToKeepGeom[gCtr] = getResult("Label",i);
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
	Table.rename(imgNameWOExt+"hoechst-intensity-measurements", "Results");
	
	// collect list of nuclei label that pass the intensity-based criteria
	labelToKeepInt = newArray(nResults);
	iCtr = 0;
	for (i=0; i<nResults; i++){
		keepCurrent = 1;
		s = getResult("Skewness",i);
		if(s > maxSkewness){
			keepCurrent = 0;
		}
		
		k = getResult("Kurtosis",i);
		if(k > maxKurtosis){
			keepCurrent = 0;
		}
		cv = getResult("StdDev",i)/getResult("Mean",i);
		if((cv<minCV) || (cv>maxCV)){
			keepCurrent = 0;
		}
		
		if(keepCurrent == 1){
			labelToKeepInt[iCtr] = getResult("Label",i);
			iCtr = iCtr+1;
		}
	}
	labelToKeepInt = Array.trim(labelToKeepInt, ctr);
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
	run("Select Label(s)", labelString);
	print("Selected "+ counter +" objects with following labels: "+LabelString);
	rename(outputImageTitle);
	
	return;
}