macro "test normalizePixelValuesToBitDepth"{
	imgName = getTitle();
	normalizePixelValuesToBitDepth(imgName);
	
}

// takes each channel of the image and applies an affine transformation 
// so that the min value of the stack (or 2D image) in each channel is zero and the max
// value is the bit depth (255 for 8-bit and 65535 for 16 bit).
// replaces the original image with the renormalized one.
function normalizePixelValuesToBitDepth(imgName){
	selectWindow(imgName);
	bd = bitDepth();
	if((bd!= 8) && (bd!=16)){
		print("img "+imgName+" is neither 8 nor 16 bit, cannot normalize.");
		return;
	}
	// set range of renormalized image to entire bitdepth
	bitMax = pow(2, bd)-1;
	
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	if(C>1){
		print("image has "+C+" channels. will normalize each channel to bit depth.");
	}
	for (i = 1; i <= C; i++) {
		selectWindow(imgName);
		run("Duplicate...", "duplicate channels="+i);
		rename("tmpChannel");
		
		// find the min and max of the image
		if(sizeZ>1){
			run("Z Project...", "projection=[Max Intensity]");
			getStatistics(area, mean, min, imgMax, std);
			close();
			run("Z Project...", "projection=[Min Intensity]");
			getStatistics(area, mean, imgMin, max, std);
			close();
		}else {
			getStatistics(area, mean, imgMin, imgMax, std);
		}
		selectWindow("tmpChannel");
		// set min pixel(s) to zero
		run("Subtract...", "value="+imgMin+" stack");
		
		// normalize to bit depth
		if(imgMin < imgMax){
			multFactor = bitMax/(imgMax-imgMin);
			if(sizeZ>1){
				run("Multiply...", "value="+multFactor+" stack");	
			}else {
				run("Multiply...", "value="+multFactor);
			}
		}
		
		if(i==1){
			selectWindow("tmpChannel");
			rename("renormImg");
		}else{
			addChannelToImg("renormImg","tmpChannel","renormImg",0);
		}
	}
	
	// replace original image with new one
	close(imgName);
	selectWindow("renormImg");
	rename(imgName);
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
