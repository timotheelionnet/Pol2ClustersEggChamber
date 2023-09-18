Prerequisite for ImageJ integration into Matlab pipeline

As per https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab

Installation
1) Put mij.jar into the java directory of Matlab (e.g for Window Machine 'C:\Program Files\MATLAB\R2009b\java\').
2) Copy also ij.jar (ImageJ) in the java directory of Matlab. Get this file from the ImageJ website: http://rsb.info.nih.gov/ij/
3) Extend the java classpath to mij.jar, e.g using the Matlab command: javaaddpath 'C:\Program Files\MATLAB\R2009b\java\mij.jar'.
4) Extend the java classpath to ij.jar, e.g using the Matlab command: javaaddpath 'C:\Program Files\MATLAB\R2009b\java\ij.jar'.
5) Start MIJ by running the Matlab command: MIJ.start; or MIJ.start("imagej-path");

Java Source Code
http://bigwww.epfl.ch/sage/soft/mij/MIJ.java


License info:
Daniel Sage, Dimiter Prodanov, Jean-Yves Tinevez and Johannes Schindelin, "MIJ: Making Interoperability Between ImageJ and Matlab Possible", ImageJ User & Developer Conference, 24-26 October 2012, Luxembourg.
http://bigwww.epfl.ch/publications/sage1205.html

Cbrewer


MorpholibJ

