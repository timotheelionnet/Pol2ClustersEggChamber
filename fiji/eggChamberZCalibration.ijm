macro "eggChamberZCalibration"{
	setBatchMode(true);
	savePath = "";
	
	inputWindowName = getTitle();
	measurementChannel = 1;
	deTrendedIntTitle = "detrended";
	eggChamberIntensityMeasurementAllChannels2(
		savePath,"eggChamber",inputWindowName,measurementChannel,deTrendedIntTitle);
	
	selectWindow(inputWindowName);
	run("Reslice [/]...", "output=0.300 start=Top");
	run("Z Project...", "projection=[Max Intensity]");

	selectWindow(deTrendedIntTitle);
	run("Reslice [/]...", "output=0.300 start=Top");
	run("Z Project...", "projection=[Max Intensity]");
	
	setBatchMode("exit and display");

}

// takes an image (inputWindowName) generates a global mask for the entire egg chamber (or the whole image)
// using the intensity in measurementChannel as a basis (and using some filtering/thresholding), 
// subtracts to each channel the minimum intensity for that channel across the entire stack.
// 
// set maskType to either "eggChamber" or "wholeImg" to decide whether the trend is based on the 
// full image or just the egg chamber
function eggChamberIntensityMeasurementAllChannels2(savePath,maskType,
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


