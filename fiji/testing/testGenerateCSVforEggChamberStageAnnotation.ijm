macro "testGenerateCSVforEggChamberStageAnnotation"{
	outFolder = "/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img";
	subDirName = "Ctrl";
	fileName = "8-Ctrl-646,MPM2-488,Ser5ph-Cy3-1zCorr";
	eggChamberSegFolderName = "eggChamberSEG";
	eggChamberSuffix = "_eggChamber";
	generateEggChamberCSV(outFolder,subDirName,
						fileName,eggChamberSegFolderName,eggChamberSuffix);
}



function generateEggChamberCSV(outFolder,subDirName,
						fileName,eggChamberSegFolderName,eggChamberSuffix){
	// make sure directory names all have a "/" at the end
	if (endsWith(outFolder, "/")!=true){
		outFolder = outFolder + "/";
	}
	if ((lengthOf(subDirName)!=0) && (endsWith(subDirName, "/")!=true)){
		subDirName = subDirName + "/";
	}
	if (endsWith(fileName, "/")!=true){
		fileNameWOExt = fileName;
		fileName = fileName + "/";
	}else{
		fileNameWOExt = substring(fileName, 0, lengthOf(fileName)-1);
	}
	if (endsWith(eggChamberSegFolderName, "/")!=true){
		eggChamberSegFolderName = eggChamberSegFolderName + "/";
	}
	
	
	ecDir = outFolder+subDirName+fileName+eggChamberSegFolderName;
	
	// find all the egg chamber files in the folder and extract their indices
	fList =  getFileList(ecDir);
	ctr = 0;
	ecIdx = newArray(lengthOf(fList));
	for(i=0;i<lengthOf(fList);i++){
		if(fList[i].contains(eggChamberSuffix)){
			extIdx = lastIndexOf(fList[i], ".");
			ecIdx[ctr] = substring(fList[i], 
						indexOf(fList[i],eggChamberSuffix)+lengthOf(eggChamberSuffix), 
						extIdx);			
			ctr = ctr+1;			
		}
	}
	ecIdx = Array.trim(ecIdx,ctr);
	
	// generate a text file in the eggChamber dir that will store the indices of
	// each egg chamber and a place holder
	csvFilePath = ecDir+"eggChamberStages.csv";
	placeHolderColumn = newArray(ctr);
	generate2columnCSVfromArrays(csvFilePath,ecIdx,placeHolderColumn,
								"eggChamberID","eggChamberStage");
	
}

// saves 2 arrays as a csv file in 2 columns.
function generate2columnCSVfromArrays(filePath,col1Array,col2Array,header1,header2){
	// Define the value of n
    
    if (lengthOf(col1Array) != lengthOf(col2Array)){
    	print("2 arrays have different sizes, cannot save as csv table");
    	return;
    }
    // Create an empty string to store the CSV data
    csvData = header1+","+header2+"\n";
    // Loop through the numbers 1 to n and add them to the CSV data with a zero in the second column
    for (i = 0; i < lengthOf(col1Array); i++) {
        csvData = csvData + col1Array[i] + ","+col2Array[i] + "\n";
    }
    
    // Save the CSV file
   f = File.open(filePath);
   print(f,csvData);
   File.close(f);
}