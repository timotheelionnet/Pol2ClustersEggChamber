macro "eggChamberZCalibration"{
	setBatchMode(true);
	savePath = "";
	inputWindowName = "Ser2Pmpm-2-f2006-1.tif";
	measurementChannel = 1;
	eggChamberIntensityMeasurementAllChannels2(savePath,inputWindowName,measurementChannel);
	setBatchMode("exit and display");
}



// takes an image (inputWindowName) generates a global mask for the entire egg chamber
// using the intensity in measurementChannel as a basis (and using some filtering/thresholding), 
// then measures intensity/neighbors etc in all channels using the new mask as an ROI
// and saves the resulting tables in the folder savePath with names C1_eggChamberInt.csv, C2_eggChamberInt.csv etc
function eggChamberIntensityMeasurementAllChannels2(savePath,inputWindowName,measurementChannel){
	selectWindow(inputWindowName);
	
	// generate a mask that encompasses the egg chamber, to be used for measurements
	selectWindow(inputWindowName);
	run("Duplicate...", "duplicate channels="+measurementChannel);
	run("Bandpass Filter...", "filter_large=2000 filter_small=100 suppress=None tolerance=5 process");
	run("Auto Threshold", "method=Default white stack");
	rename("eggChamberMask");
	
	// computes the median across the egg Chamber mask as a function of z and c
	// and the minimum value across x,y as a function of z and c
	// outputs are two 16-bit 2D images (nz nc)
	//computeMedianAndMinZTrends(inputWindowName,"eggChamberMask","medianEggChamber","minWholeImg");
	computeMeanAndMinZTrends(inputWindowName,"eggChamberMask","meanEggChamber","minWholeImg");
	
	// compute in each channel the linear trend of the eggchamber median Z intensity 
	// (after subtracting the min intensity at each slice);
	subtractMinFromTrend("meanEggChamber","minWholeImg","corrMean");
	computeLinearTrendAlongX("corrMean","linTrend");

	// correct the intensity at each slice by the linear trend
	subtractBackgroundSliceWise(inputWindowName,"minWholeImg","bgCorrInt");
	correctIntensityForZTrend("bgCorrInt","linTrend","deTrendedInt");
}

//
function subtractMinFromTrend(trendImgTitle,minImgTitle,outputImgTitle){
	selectWindow(trendImgTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	// find the minimum of each row
	minArray = newArray(sizeY);
	selectWindow(minImgTitle);
	for (j = 0; j < sizeY; j++) {
		for (i = 0; i < sizeX; i++) {
			if(i==0){
				minArray[j] = getPixel(i,j);
			}else{
				if(getPixel(i,j)<minArray[j]){
					minArray[j] = getPixel(i,j);
				}
			}
		}
	}
	
	//duplicate trendImgTitle
	selectWindow(trendImgTitle);
	run("Select All");
	run("Duplicate...", "duplicate");
	rename(outputImgTitle);
	for (j = 0; j < sizeY; j++) {
		for (i = 0; i < sizeX; i++) {
			setPixel(i, j, getPixel(i, j) - minArray[j]);
		}
	}
}

// from an input hyperstack inputImg (nx ny nz nc), and a 2D image minTrendTitle (nz,nc) holding the min intensity vs z by channel 
// generates a corrected hyperstack outputImg (nx ny nz nc) where the intensity in each slice is corrected
// Icorr(x,y,z,c) = Iinput(x,y,z,c) * linTrend(0,c) / linTrend(z,c)
function subtractBackgroundSliceWise(inputImgTitle,minTrendTitle,outputImgTitle){
	selectWindow(inputImgTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	//duplicte input hyperstack
	run("Select All");
	run("Duplicate...", "duplicate");
	rename(outputImgTitle);
	run("32-bit");

	// subtract background at each slice
	curTrend = newArray(C);
	for(izs = 0; izs<sizeZ; izs++) {
		for(ic = 0; ic<C; ic++) {
			selectWindow(minTrendTitle);
			curMin = getPixel(izs, ic);
			
			selectWindow(outputImgTitle);
			Stack.setPosition(ic+1, izs+1, 1);
			run("Subtract...", "value="+ curMin +" slice");
		}
	}
}	


// from an input hyperstack inputImgTitle (nx ny nz nc), and a 2D image linTrendTitle (nz,nc) holding the linear trend vs z by channel 
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

// from hyperstack hs (nx ny nz nc), and a mask (nx nz nz 1 8-bit stack) at each z-slice
// computes the median across the mask as a function of z and c
// and the minimum value across x,y as a function of z and c
// outputs are two 16-bit 2D images (nz nc), one called medianTrendTitle (holdig the median)
// one called minTrendTitle (holding the min).
function computeMedianAndMinZTrends(hsTitle,maskTitle,medianTrendTitle,minTrendTitle){
	selectWindow(hsTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	newImage(medianTrendTitle, "16-bit", sizeZ, C, 1);
	newImage(minTrendTitle, "16-bit", sizeZ, C, 1);
	
	meanChannels = newArray
	for(izs = 0; izs<sizeZ; izs++) {
		roiManager("reset");
	    selectWindow(maskTitle);
	    Stack.setPosition(1, izs+1, 1);
	    run("Create Selection");
	    roiManager("Add");
		for (i = 0; i < C; i++) {
		    selectWindow(hsTitle);
		    Stack.setPosition(i+1, izs+1, 1);
		    run("Select All");
		    getStatistics(area, mean, minImg, max, std, histogram);
			selectWindow(minTrendTitle);
			setPixel(izs, i, minImg);
			
			selectWindow(hsTitle);
			roiManager("Select", 0);
			medianEC = getValue("Median");
		    selectWindow(medianTrendTitle);
			setPixel(izs, i, medianEC);
		}
	}
}

// from hyperstack hs (nx ny nz nc), and a mask (nx nz nz 1 8-bit stack) at each z-slice
// computes the mean across the mask as a function of z and c
// and the minimum value across x,y as a function of z and c
// outputs are two 16-bit 2D images (nz nc), one called meanTrendTitle (holdig the mean)
// one called minTrendTitle (holding the min).
function computeMeanAndMinZTrends(hsTitle,maskTitle,meanTrendTitle,minTrendTitle){
	selectWindow(hsTitle);
	getDimensions(sizeX, sizeY, C, sizeZ, F);

	newImage(meanTrendTitle, "16-bit", sizeZ, C, 1);
	newImage(minTrendTitle, "16-bit", sizeZ, C, 1);
	
	meanChannels = newArray
	for(izs = 0; izs<sizeZ; izs++) {
		roiManager("reset");
	    selectWindow(maskTitle);
	    Stack.setPosition(1, izs+1, 1);
	    run("Create Selection");
	    roiManager("Add");
		for (i = 0; i < C; i++) {
		    selectWindow(hsTitle);
		    Stack.setPosition(i+1, izs+1, 1);
		    run("Select All");
		    getStatistics(area, mean, minImg, max, std, histogram);
			selectWindow(minTrendTitle);
			setPixel(izs, i, minImg);
			
			selectWindow(hsTitle);
			roiManager("Select", 0);
			getStatistics(area, meanEC, minImg, max, std, histogram);
		    selectWindow(meanTrendTitle);
			setPixel(izs, i, meanEC);
		}
	}
}

// takes input 2D image titled inputImg and generates a 2D image titled linTrend
// with same size where the value of each pixel (x,y) is equal the linear interpolation of its row
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


