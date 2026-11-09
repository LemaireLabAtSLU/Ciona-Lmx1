//// This script requires the plugin 3Dsuite https://mcib3d.frama.io/3d-suite-imagej/#installation

//setBatchMode(true);
dir1 = getDirectory("Choose Source Directory for image to binarise"); //source image
index = lastIndexOf(dir1,"\\");
print(dir1);
print(index);
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_Lmx-binary/";
File.makeDirectory(dir2);

dir4 = substring(dir1,0,index);
j = getFileList(dir1);
print(j.length);
for (i=0; i<j.length; i++)
{
	open(dir1 + j[i]);
	title1=getTitle();
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
	selectWindow(title1);
	getDimensions(width, height, channels, slices, frames);
	run("Split Channels");
	
	
	close("C1-"+title2+".tif");
	close("C2-"+title2+".tif");
	close("C4-"+title2+".tif");
	close("C5-"+title2+".tif");
	
	selectWindow("C3-"+title2+".tif");
	rename("Lmx1");
	run("Duplicate...", "duplicate");
	rename("Lmx1 background");

	run("Gaussian Blur...", "sigma=90 stack");
	run("Multiply...", "value=1.15 stack");

	imageCalculator("Subtract create stack", "Lmx1","Lmx1 background");
	
	close("Lmx1");
	close("Lmx1 background");
	selectWindow("Result of Lmx1");
	
	run("Kuwahara Filter", "sampling=7 stack");

	run("3D Iterative Thresholding", "min_vol_pix=200 max_vol_pix=3000 min_threshold=80 min_contrast=.1 criteria_method=VOLUME threshold_method=STEP segment_results=All value_method=10");
	
	if (isOpen("draw")) 
		{
			
		selectWindow("draw");
	
		rename("Lmx1-draw");
		
		getDimensions(width, height, channels, slices, frames);
		//print(width, height, channels, slices, frames);
		if (channels == 1) {

	
		logString = getInfo( "log" ); //Obtaining the number of Lmx Nuclei
		index2 = lastIndexOf(logString, ":");
		index3 = lengthOf(logString);
		lmxNucleiNumber = substring(logString,index2+2,index3);
		print("Number of nuclei: "+lmxNucleiNumber);	

		setThreshold(1, 65535, "raw");
		run("Make Binary", "method=Default background=Default create");

		selectWindow("MASK_Lmx1-draw");
		rename("MASK_Lmx1");
		run("3D Maxima Finder", "minimmum=0 radiusxy=1.50 radiusz=1.50 noise=100");
		selectWindow("peaks_MASK_Lmx1");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		
		selectWindow("MASK_Lmx1");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");

		run("3D Watershed Split", "binary=MASK_Lmx1 seeds=peaks_MASK_Lmx1 radius=2");
		close("EDT");
		selectWindow("Split");
		setThreshold(1, 65535, "raw");
		run("Make Binary", "method=Default background=Default create");
		close("Split");

		selectWindow("MASK_Split");
		saveAs("tiff", dir2+"Lmx1-Binary_"+title2+".tif");
		
		}
		else 
		{
		selectWindow("Lmx1-draw");
		run("Make Composite", "display=Composite");
    	run("Split Channels");
		
		for (a = 1; a <= channels; a++) 
		{
    	selectWindow("C"+a+"-Lmx1-draw");
    	setThreshold(1, 65535, "raw");
		run("Make Binary", "method=Default background=Default create");
		selectWindow("MASK_C"+a+"-Lmx1-draw");
		rename("C"+a+"MASK_Lmx1");
		
		run("3D Maxima Finder", "minimmum=0 radiusxy=1.50 radiusz=1.50 noise=100"); //Find maxima for watershed after dilation
		selectWindow("peaks_C"+a+"MASK_Lmx1");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		
		selectWindow("C"+a+"MASK_Lmx1");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");

		run("3D Watershed Split", "binary=C"+a+"MASK_Lmx1 seeds=peaks_C"+a+"MASK_Lmx1 radius=2");
		close("peaks_C"+a+"MASK_Lmx1");
		close("C"+a+"MASK_Lmx1");
		close("EDT");
		selectWindow("Split");
		setThreshold(1, 65535, "raw");
		run("Make Binary", "method=Default background=Default create");
		close("Split");

		selectWindow("MASK_Split");
		rename("C"+a+"-MASK_Lmx1_w");

		}
		args = "c1=C1-MASK_Lmx1_w";
		
			for (c = 2; c <= channels; c++) 
		{
    	args = args+" c"+c+"=C"+c+"-MASK_Lmx1_w";
    	print(args);
		}
		args = args+" create";
		run("Merge Channels...", args);
		saveAs("tiff", dir2+"Lmx1-Binary_"+title2+".tif");
		}
		} 
		else 
		{
    	print("Skipping " + title1 + ": Lmx1 is not detected.");
		}


	close("Lmx1-Binary_"+title2+".tif");
	close("Result of Lmx1");

	run("Close All");
}
selectWindow("Log");
saveAs("Text", dir4+"_BinaryLogFile.csv");
//setBatchMode(false);

print("\nFinish");