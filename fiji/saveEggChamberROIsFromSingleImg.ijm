macro "saveEggChamberROIsFromSingleImage"{
	outFolder = getDirectory("choose the output directory");
	originalImgTitle = getTitle();	
	getDimensions(sizeX, sizeY, C, sizeZ, F);
	
	k = lastIndexOf(originalImgTitle, ".");
	
	n = roiManager("count");
	for (i = 0; i < n; i++) {
		selectWindow(originalImgTitle);
		roiManager("Select", i);
		run("Create Mask");
		run("16-bit");
		run("Divide...", "value=255");
		idx = i+1; // make sure that indices start at 1, not 0
		run("Multiply...", "value="+idx);
		setMinAndMax(0, idx);
		imgName = substring(originalImgTitle,0,k-1)+"_eggChamber"+idx+".tif";
		rename(imgName);
		save(outFolder+imgName);
		close();
	}
}
