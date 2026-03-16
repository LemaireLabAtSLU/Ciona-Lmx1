dir1 = getDirectory("Image Folder");
dir2 = getDirectory("Binary Folder");
index = lastIndexOf(dir1,"\\");
tempdir = substring(dir1,0,index);
dir3 = tempdir+"_Stat/";
File.makeDirectory(dir3);
j = getFileList(dir1);
for (i=0; i<j.length; i++)
{
	open(dir1 + j[i]);
	title1=getTitle();
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
	getDimensions(w, h, channels, slices, frames); 
	print("Chanels: "+channels);
	print("Slices: "+slices);
	print("Frames: "+frames);
	run("Specify...", "width=40 height=40 x=256 y=256 slice=1");
 	waitForUser("select background area");
 	run("Duplicate...", "duplicate channels=2");
 	title3=getTitle();
	run("Z Project...", "projection=[Sum Slices]");
	run("Measure");
	totalBackground = getResult("Mean", 0);
	background = totalBackground/slices;
	print("Background: "+background);
	close(title3);
	run("Clear Results");
	selectWindow(title1);	
	run("Select All");
	run("Duplicate...", "duplicate channels=2");
	rename("C2-"+title2);
	run("Subtract...", "value="+background+" stack");
	open(dir2 + "LmxBinary_"+title1);
	title4 = getTitle();
	run("3D Objects Counter", "threshold=128 slice=15 min.=100 max.=5000 objects statistics summary");
	selectWindow("Statistics for "+title4);
	saveAs("Results", dir3+"Nucleus_"+title2+".csv");
	close("Nucleus_"+title2+".csv");
	close("Objects map of LmxBinary_"+title1);
	selectWindow("LmxBinary_"+title1);
	run("Dilate (3D)", "iso=255");	
	run("Dilate", "stack");
	run("Dilate", "stack");
	totalIntensity = 0;
for (b=0; b<slices; b++)
{
selectWindow("LmxBinary_"+title1);
Stack.setPosition(1, b+1, 1);
run("Create Selection");
run("Measure");	
area = getResult("Area", b);
//print("area: "+area);
if (area == 0.00)
{
sliceIntensity = 0;
}
else
{
Table.deleteRows(b,b);	
selectWindow("C2-"+title2);
Stack.setPosition(1, b+1, 1);
run("Restore Selection");
run("Measure");	
sliceIntensity = getResult("IntDen",b);
totalIntensity = totalIntensity+ sliceIntensity;
}
}
print("Total intensity: "+totalIntensity);
selectWindow("Results");
saveAs("Results", dir3+"Intensity_"+title2+".csv");
run("Clear Results");
run("Close All");
}


