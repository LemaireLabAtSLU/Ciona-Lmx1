//setBatchMode(true);
dir1 = getDirectory("Choose Source Directory for image to binarise"); //source image
index = lastIndexOf(dir1,"\\");
print(index);
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_mask8/";
File.makeDirectory(dir2);
dir3 = tempdir+"_nucleusStat8/";
File.makeDirectory(dir3);
j = getFileList(dir1);
for (i=0; i<j.length; i++)
{
	open(dir1 + j[i]);
	title1=getTitle();
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
	selectWindow(title1);
	run("Split Channels");
	close("C1-"+title2+".tif");
	close("C2-"+title2+".tif");
	close("C4-"+title2+".tif");
	close("C5-"+title2+".tif");
	selectWindow("C3-"+title2+".tif");
	run("Duplicate...", "duplicate");
	rename("background");
	getDimensions(width, height, channels, slices, frames);
//	for (b = 0; b < slices; b++)
//	{
//    	Stack.setPosition(1, b+1, 1);
//    	setThreshold(10000, 65535, "raw");
//		run("Create Selection");
//		run("Set...", "value=10000");
//    }
//    resetThreshold();
//   run("Select All");	
	run("Median...", "radius=35 stack");
//	run("Multiply...", "value=2 stack");
	selectWindow("C3-"+title2+".tif");
	rename("embryo");
//	for(j=0; j<slices; j++)
//	{
//	Stack.setPosition(1, j+1, 1);
//	run("Enhance Local Contrast (CLAHE)", "blocksize=100 histogram=256 maximum=2 mask=*None*");
//	print(j);
//	}
	imageCalculator("Subtract create stack", "embryo","background");
	close("embryo");
	close("background");
	selectWindow("Result of embryo");
	run("Smooth", "stack");
	run("Smooth", "stack");
	run("Smooth", "stack");
	//run("Smooth", "stack");
	//title3 = getTitle();
	//print(title3);
	//selectWindow(title3);
	setMinAndMax(0, 30000); 
	setOption("ScaleConversions", true);
	run("8-bit");
	run("Auto Local Threshold", "method=Bernsen radius=15 parameter_1=0 parameter_2=0 white stack");
	run("Make Binary", "method=Default background=Default create");
	close("Result of embryo");
	selectWindow("MASK_Result of embryo");
	run("Fill Holes", "stack");
	run("Erode", "stack");
	run("Dilate", "stack");
	run("Dilate", "stack");
	run("Erode", "stack");
	run("Watershed", "stack");
	run("Options...", "iterations=2 count=3 do=Open stack");
	run("Options...", "iterations=3 count=3 do=Close stack");
//	run("Options...", "iterations=3 count=3 do=Dilate stack");
	run("Watershed", "stack");
	title4 = getTitle();
	print(title4);
	run("3D Objects Counter", "threshold=128 slice=15 min.=100 max.=5000 objects statistics summary");
	selectWindow("Statistics for "+title4);
	saveAs("Results", dir3+"NucleusResults_"+title2+".csv");
	close("NucleusResults_"+title2+".csv");
	selectWindow(title4);
	saveAs("tiff", dir2+"LmxBinary_"+title2+".tif");	
	run("Close All");
}
selectWindow("Log");
saveAs("Results", dir3+"LogFile.csv");
//setBatchMode(false);