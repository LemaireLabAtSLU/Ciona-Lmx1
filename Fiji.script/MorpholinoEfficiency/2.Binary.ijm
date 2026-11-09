//// This script requires the plugin 3Dsuite https://mcib3d.frama.io/3d-suite-imagej/#installation


setBatchMode(true);
dir1 = getDirectory("Choose Source Directory for image to binarise"); //source image
index = lastIndexOf(dir1,"\\");
print(dir1);
print(index);
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_Lmx-binary/";
File.makeDirectory(dir2);
dir3 = tempdir+"_Msx-binary/";
File.makeDirectory(dir3);
dir4 = substring(dir1,0,index);
j = getFileList(dir1);
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
	
	
	close("C3-"+title2+".tif");
	
	selectWindow("C1-"+title2+".tif");
	rename("Lmx1");
	run("Duplicate...", "duplicate");
	rename("Lmx1 background");

	run("Gaussian Blur...", "sigma=35 stack");
	run("Multiply...", "value=1.25 stack");

	imageCalculator("Subtract create stack", "Lmx1","Lmx1 background");
	close("Lmx1");
	close("Lmx1 background");
	selectWindow("Result of Lmx1");
	run("Kuwahara Filter", "sampling=7 stack");
	
	run("3D Iterative Thresholding", "min_vol_pix=100 max_vol_pix=1000 min_threshold=25 min_contrast=0 criteria_method=ELONGATION threshold_method=STEP segment_results=All value_method=10");

	
		if (isOpen("draw")) 
		{
		selectWindow("draw");
	
		rename("Lmx1-draw");
	
		logString = getInfo( "log" ); //Obtaining the number of Lmx Nuclei
		index2 = lastIndexOf(logString, ":");
		index3 = lengthOf(logString);
		lmxNucleiNumber = substring(logString,index2+2,index3);
		print("Number of nuclei: "+lmxNucleiNumber);	

		setThreshold(1, 65535, "raw");
		run("Make Binary", "method=Default background=Default create");

		selectWindow("MASK_Lmx1-draw");
		rename("MASK_Lmx1");

		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");
		run("Dilate (3D)", "iso=255");

	
		title4 = getTitle();
		print(title4);
		saveAs("tiff", dir2+"Lmx1-Binary_"+title2+".tif");
		} 
		else 
		{
    	print("Skipping " + title1 + ": it does not express Lmx1.");
		}

	close("Lmx1-Binary_"+title2+".tif");
	close("Result of Lmx1");

	selectWindow("C2-"+title2+".tif");
	rename("Msx");
	run("Duplicate...", "duplicate");
	rename("Msx background");
	run("Gaussian Blur...", "sigma=85 stack");
	imageCalculator("Subtract create stack", "Msx","Msx background");
	close("Msx");
	close("Msx background");
	selectWindow("Result of Msx");
	run("Kuwahara Filter", "sampling=7 stack");
	run("Auto Threshold", "method=Huang2 white show stack use_stack_histogram");
	run("Convert to Mask", "background=Dark list create");
	rename("Mask-Msx");
	run("Fill Holes", "stack");
	run("Options...", "iterations=3 count=2 do=Open stack");
	run("Dilate (3D)", "iso=255");
	run("Dilate (3D)", "iso=255");
	run("3D Binary Close Labels", "radiusxy=6 radiusz=3 operation=Close");
	run("Convert to Mask", "background=Dark list create");
	saveAs("tiff", dir3+"Msx-Binary_"+title2+".tif");

	run("Close All");
}
selectWindow("Log");
saveAs("Text", dir4+"_BinaryLogFile.csv");
setBatchMode(false);

print("\nFinish");