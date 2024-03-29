macro "batch detrend"{
	
	//root folder where all data is
	input = getDirectory("choose the input directory");
	
	//where the analysis images will go 
	output = getDirectory("choose the output directory");

	if (endsWith(output,File.separator) == false ){
		output = output + File.separator;
	}

	if (endsWith(input,File.separator) == false ){
		input = input + File.separator;
	}
			
	//batch mode so nothing is displayed to accelerate execution
	setBatchMode(true); 
	print("input dir: "+input);
	print("output dir: "+output);
	print("de-trending intensity along Z ...");

	//create analysis directories if they do not exist
	
	//collect file and folder names from root dir
	tmplist = getFileList(input);
	
	//declare arbitrarily large file list
	list = newArray(10000);
	
	//counter of the image files found in the directory and subdirectories
	ctr = 0;
	
	//search through arborescence for tif images
	for (k=0; k<tmplist.length; k++) {
	 	
	    if (endsWith(tmplist[k], "/") || endsWith(tmplist[k], "\\") ){
			tmplist2 = getFileList(input+tmplist[k]);
			print(tmplist[k]);	
			for (l=0; l<tmplist2.length; l++) {
				print(tmplist2[l]);	
				if (endsWith(tmplist2[l], "/") || endsWith(tmplist2[l], "\\") ){
					tmplist3 = getFileList(input+tmplist[k]+tmplist2[l]);
					for (m=0; m<tmplist3.length; m++) {	
						if( (indexOf(tmplist3[m], ".tif") != -1) || (indexOf(tmplist3[m], ".lsm") != -1) || (indexOf(tmplist3[m], ".nd2") != -1) || (indexOf(tmplist3[m], ".czi") != -1)){
		        			list[ctr] = input+tmplist[k]+tmplist2[l]+tmplist3[m];
		        			//print(list[ctr]);
		        			ctr++;
	        			}
					}
				}
				if((indexOf(tmplist2[l], ".tif") != -1) || (indexOf(tmplist2[l], ".lsm") != -1) || (indexOf(tmplist2[l], ".nd2") != -1) || (indexOf(tmplist2[l], ".czi") != -1)){
		        	list[ctr] = input+tmplist[k]+tmplist2[l];
		        	//print(list[ctr]);
		        	ctr++;
	        	}
			}	
		}  
	
	    if((indexOf(tmplist[k], ".tif") != -1) || (indexOf(tmplist[k], ".lsm") != -1) || (indexOf(tmplist[k], ".nd2") != -1) || (indexOf(tmplist[k], ".czi") != -1)){
		      list[ctr] = input+tmplist[k];
		      //print(list[ctr]);
		      ctr++;
	    }
	         
	 }
	 list = Array.trim(list,ctr);
	 
	 // fixing bug seen on a PC that some subfolders have their last file separator a slash instead of a backslash
	 // so force changing all slashes or backslashes to the right file separator
	 switchFileSep = 1;
	 if(switchFileSep == 1){
	 	for (k=0; k<list.length; k++) {
	 		if (File.separator == "\\"){
	 			list[k] = replace(list[k], "/", File.separator);
	 		}else {
	 			list[k] = replace(list[k], "\\", File.separator);
	 		}
	 	}
	 }
	 
	 print("list of files to treat:");
	 for (k=0; k<list.length; k++) {
	 	print(list[k]);
	 }
	 
	print("starting file processing...");
	//now that file names for the conditions are all combined, lets do the analysis
	for (i = 0; i < list.length; i++){
			
        if ((indexOf(list[i], ".tif") != -1) || (indexOf(list[i], ".lsm") != -1) || (indexOf(list[i], ".nd2") != -1) || (indexOf(list[i], ".czi") != -1)) { // making sure it is an image
	        //open image, split channels and asve each
	        open(list[i]); 
        	originalImgTitle = getTitle();	
        	print("processing file "+ originalImgTitle);
			n = indexOf(originalImgTitle,".");
			strname = substring(originalImgTitle,0,n);
			zcorrImgTitle = strname+"zCorr.tif";
			eggChamberIntensityMeasurementAllChannels2("eggChamber",
				originalImgTitle,1,zcorrImgTitle);
			
			// remove the input folder from the original image file
			n1 = indexOf(list[i], input);
			if (n1 != -1){
				
			}else{
				print("File ",list[i]," does not seem to belong to input folder");
			}
			n1 = indexOf(list[i], input);
			if (n1 != -1){
				outSubfolder = substring(list[i],n1+ lengthOf(input));
				print("outSubfolder: "+outSubfolder);
			}else{
				print("File ",list[i]," does not seem to belong to input folder");
			}
			// add any subfolder in front of image name if needed.
			n2 = lastIndexOf(outSubfolder,File.separator);
			if (n2!=-1){
				outSubfolder = substring(outSubfolder, 0,n2+1);
				print("outSubfolder: "+outSubfolder);
			}	
			print("saving to "+ output + outSubfolder+ zcorrImgTitle +" ...");

			// making sure directory exists or create it
			subDirList = split(outSubfolder, File.separator);
			curDirToCheck = output;
			for (j = 0; j < subDirList.length; j++) {
				if(endsWith(curDirToCheck, File.separator)== false){
					curDirToCheck = curDirToCheck + File.separator;
				}
				curDirToCheck = curDirToCheck+subDirList[j];
				if(File.exists(curDirToCheck)){
					print("subDir "+curDirToCheck+" exists");
				}else {
					print("subDir "+curDirToCheck+" doesn't exist, creating it...");
					File.makeDirectory(curDirToCheck);
				}
			}
			
			saveAs("Tiff", output + outSubfolder + zcorrImgTitle);
			close(originalImgTitle);
			close(zcorrImgTitle);
			close(zcorrImgTitle+".tif"); // I think saving the image appends .tif to its name. 
			run("Collect Garbage");
        }
	}
	print("done.");
	setBatchMode("exit and display");		
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