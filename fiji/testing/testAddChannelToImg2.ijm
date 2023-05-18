macro "testAddChannelToImg2"{
	run("Close All");
	open("/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/Ctrl/8-Ctrl-646,MPM2-488,Ser5ph-Cy3-1001zCorr/nucTIF/nuc1_masks_plasmCorr.tif");
	imgSource = "nuc1_masks_plasmCorr.tif";
	
	// uncomment this next block of code to test single color images
	/*
	selectWindow("nuc1_masks_plasmCorr.tif");
	run("Duplicate...", "title=c2Duplicate duplicate channels=2");
	imgSource = "c2Duplicate";
	*/
	//
	
	selectWindow("nuc1_masks_plasmCorr.tif");
	run("Duplicate...", "title=c1Duplicate duplicate channels=1");
	channelSource = "c1Duplicate";
	
	newImgName = "res";
	keepSourceImgs = 1;
	addChannelToImg2(imgSource,channelSource,newImgName,keepSourceImgs);
}


// appends an extra color channel (channelSource) to a hyperstack (imgSource) and renames the result newImgName 
function addChannelToImg2(imgSource,channelSource,newImgName,keepSourceImgs){

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
	
	if(keepSourceImgs == 1){
		keepString = " keep";
	}else{
		keepString = "";
	}
	
	strResStart = "  title=" + newImgName + keepString + " open ";
	strRecomposeStart = "  title=" + imgSource + " open ";
	
	selectWindow(imgSource);
	getDimensions(w, h, c, z, f);
	if(c==1){
		argumentString1 = "image1="+imgSource+" image2="+channelSource+" image3=[-- None --]";
	}else{
		run("Split Channels");
		argumentString0 = "";
		for(i=1;i<=c;i++){
			argumentString0 = argumentString0 + "image"+i+"=C"+i+"-"+imgSource+" ";
		}
		i = c+1;
		argumentString1 = argumentString0 + "image"+i+"="+channelSource;
		argumentString0 = argumentString0 +" image"+i+"=[-- None --]";
		i = i+1;
		argumentString1 = argumentString1 +" image"+i+"=[-- None --]";
	}

	if(keepSourceImgs == 1){
		// add channel to split channels from original image, keep sources
		run("Concatenate...", strResStart + argumentString1);
		print("s1: "+ strResStart + argumentString1);	
		rename(newImgName);
		run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
		if(c+1>=8){
			Stack.setDisplayMode("color");
		}
		if(c!=1){
			//recombine split channels into original image
			print("s0: "+ strRecomposeStart + argumentString0);
			run("Concatenate...", strRecomposeStart + argumentString0);
			run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
			if(c>=8){
			Stack.setDisplayMode("color");
		}
			for(i=1;i<=c;i++){
				close("C"+i+"-"+imgSource);
			}
		}
	}else{
		print("s1: "+ strResStart + argumentString1);	
		run("Concatenate...", strResStart + argumentString1);
		run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
		if(c+1>=8){
			Stack.setDisplayMode("color");
		}
		rename(newImgName);
	}
}