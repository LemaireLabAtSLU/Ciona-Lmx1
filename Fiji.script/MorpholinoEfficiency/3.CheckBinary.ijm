dir1 = getDirectory("Image Folder");

index = lastIndexOf(dir1,"\\");
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_Lmx-binary/"; //directory of the Lmx binary files
dir3 = tempdir+"_Msx-binary/"; //directory of the Msx binary files

dir4 = tempdir+"_LmxCheckedBinary/";
File.makeDirectory(dir4);
dir5 = tempdir+"_MsxCheckedBinary/";
File.makeDirectory(dir5);

j = getFileList(dir1);
for (i=0; i<j.length; i++)
{
	open(dir1 + j[i]);
	title1=getTitle();
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
	//Check Lmx binary
	if (File.exists(dir2 + "Lmx1-Binary_"+title1)) //Verify that Lmx binary exist
	{
	open(dir2 + "Lmx1-Binary_"+title1);
	getDimensions(width, height, channels, slices, frames);
	if (channels == 1) 
	{
	selectWindow("Lmx1-Binary_"+title1);
	rename("LmxBinary");
	run("16-bit");
	run("Invert LUT");
	run("Divide...", "value=2 stack");
	selectWindow(title1);
	run("Split Channels");
	close("C3-"+title1);
	selectWindow("C2-"+title1);
	rename("Msx");
	selectWindow("C1-"+title1);
	rename("Lmx");
	run("Merge Channels...", "c1=LmxBinary c2=Lmx create keep");
	Stack.setChannel(1);
	setMinAndMax(0, 500);
	Stack.setChannel(2);
	setMinAndMax(0, 800);
	waitForUser("Check Binary");
	run("Split Channels");
	selectWindow("C1-Composite");
	setThreshold(128, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir4+"LmxBinary_"+title2+".tif");
	run("Close All");
	}
	else 
	{ 
		//the Lmx binary has more than one channel
	selectWindow(title1);
	run("Split Channels");
	close("C3-"+title1);
	close("C2-"+title1);
	selectWindow("C1-"+title1);
	rename("Lmx");
	
	selectWindow("Lmx1-Binary_"+title1);
	rename("LmxBinary");
	run("16-bit");
	run("Divide...", "value=2 stack");
	
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+slices+" frames=1 display=Color");
	
	run("Split Channels");
	
	
	run("Concatenate...", "all_open open");
	run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	rename("LmxBinary");
	Stack.setChannel(1);
	setMinAndMax(0, 500);
	
	for (k = 0; k < channels; k++) 
	{
	Stack.setChannel(k+2);
	setMinAndMax(0, 250);
	}
	
	title = "Check and Select binary";
  	msg = "Remove non specific segmentation.";
	waitForUser(title, msg);
	Dialog.create("Specify Channel to duplicate");
	Dialog.addNumber("Selected channel for Lmx binary", 1);
	Dialog.show();

	channel_num = Dialog.getNumber();
	run("Duplicate...", "duplicate channels="+channel_num);
	
	setThreshold(128, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir4+"LmxBinary_"+title2+".tif");
	run("Close All");
	}
	} 
	else
	{ 
		//Create a binary file for Lmx
		print("No Lmx binary found");
		selectWindow(title1);
		getDimensions(width, height, channels, slices, frames);
		newImage("HyperStack", "16-bit grayscale-mode", 1024, 1024, 1, slices, 1);
		setThreshold(128, 65535, "raw");
		run("Make Binary", "method=Default background=Default create");
		saveAs("tiff", dir4+"LmxBinary_"+title2+".tif");
		run("Close All");
			

	}
	//Check Msx binary
open(dir3 + "Msx-Binary_"+title1);
	open(dir1 + title1);
	run("Split Channels");
	close("C3-"+title1);
	close("C1-"+title1);
	selectWindow("C2-"+title1);
	rename("Msx");
	selectWindow("Msx-Binary_"+title1);
	rename("MsxBinary");
	run("16-bit");
	run("Invert LUT");
	run("Divide...", "value=2 stack");

	run("Merge Channels...", "c1=MsxBinary c2=Msx create keep");
	Stack.setChannel(1);
	setMinAndMax(0, 500);
	Stack.setChannel(2);
	setMinAndMax(0, 800);
	waitForUser("Check Binary");
	run("Split Channels");
	selectWindow("C1-Composite");
	setThreshold(128, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir5+"MsxBinary_"+title2+".tif");
	run("Close All");
}
	