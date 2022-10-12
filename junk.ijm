macro ""{
	setBatchMode(true);
	imageCalculator("Subtract create 32-bit", "medianEggChamber","minWholeImg");
	rename("corrMedian");
	computeLinearTrendAlongX("corrMedian","linTrend");
	setBatchMode("exit and display");

}

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