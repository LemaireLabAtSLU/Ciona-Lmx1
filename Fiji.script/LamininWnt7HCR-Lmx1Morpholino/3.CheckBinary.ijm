dir1 = getDirectory("Image Folder");

index = lastIndexOf(dir1,"\\");
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_Lmx-binary/"; //directory of the Lmx binary files


dir4 = tempdir+"_Lmx_blin/"; //folder for binary of Lmx+ nuclei of blineage
File.makeDirectory(dir4);
dir5 = tempdir+"_Lmx_alin/"; //folder for binary of Lmx+ nuclei of alineage
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
	close("C1-"+title1);
	close("C2-"+title1);
	close("C4-"+title1);
	close("C5-"+title1);
	
	selectWindow("C3-"+title1);
	rename("Lmx");
	run("Merge Channels...", "c1=LmxBinary c2=Lmx create keep");
	Stack.setChannel(1);
	setMinAndMax(0, 500);
	Stack.setChannel(2);
	setMinAndMax(0, 800);
	
	//// Cleanup Binary if necessary
	waitForUser("Check Binary,\nSplit nuclei if necessary or\nDraw missing a or b Lmx+ dorsal neural cells");
	run("Duplicate...", "duplicate");
	selectWindow("Composite");
	rename("b-lineage");
	selectWindow("Composite-1");
	rename("a-lineage");
	
	
	/// Select b lineage nuclei
	selectWindow("b-lineage");	
	
	setTool("freehand");
	waitForUser("With FreeHand selection tool,\nselect b-lineage Lmx+ cells");
	
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
	run("Split Channels");
	close("C2-b-lineage");
	selectWindow("C1-b-lineage");
	setThreshold(1, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir4+"Lmx-b_"+title2+".tif");
	close("Lmx-b_"+title2+".tif");
	
	/// Select a lineage nuclei	
	selectWindow("a-lineage");
	setTool("freehand");
	waitForUser("With FreeHand selection tool,\nselect a-lineage Lmx+ cells");
	
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
	
	run("Split Channels");
	close("C2-a-lineage");
	selectWindow("C1-a-lineage");
	setThreshold(1, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir5+"Lmx-a_"+title2+".tif");
	close("Lmx-a_"+title2+".tif");
	
	
	run("Close All");
	}
	else 
	{ 
		//the Lmx binary has more than one channel, select the best binary
	selectWindow(title1);
	run("Split Channels");
	close("C1-"+title1);
	close("C2-"+title1);
	close("C4-"+title1);
	close("C5-"+title1);
	selectWindow("C3-"+title1);
	rename("Lmx");
	//run("Duplicate...", "duplicate");
	selectWindow("Lmx1-Binary_"+title1);
	rename("LmxBinary");
	run("16-bit");
	run("Divide...", "value=2 stack");
	
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+slices+" frames=1 display=Color");
	
	run("Split Channels");
	
	
	run("Concatenate...", "all_open open");
	run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	rename("Composite");
	Stack.setChannel(1);
	setMinAndMax(0, 500);
	
	for (k = 0; k < channels; k++) 
	{
	Stack.setChannel(k+2);
	setMinAndMax(0, 250);
	}
	run("Duplicate...", "duplicate channels=1");
	rename("Lmx-1");
	
	title = "Check and Select binary";
  	msg = "Remove non specific segmentation.";
	waitForUser(title, msg);
	run("Select All");
	print("Done");
	
	Dialog.create("Specify Channel to duplicate");
	Dialog.addNumber("Selected channel for Lmx binary", 1);
	Dialog.show();

	channel_num = Dialog.getNumber();
	run("Duplicate...", "duplicate channels="+channel_num);
	rename("LmxBinary");
	close("Composite"); 

	run("Merge Channels...", "c1=LmxBinary c2=Lmx-1 create keep");
	run("Duplicate...", "duplicate");
	selectWindow("Composite");
	rename("b-lineage");
	selectWindow("Composite-1");
	rename("a-lineage");
	
	/// Select b lineage nuclei
	selectWindow("b-lineage");	
	setTool("freehand");
	waitForUser("With FreeHand selection tool,\nselect b-lineage Lmx+ cells");
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
	run("Split Channels");
	close("C2-b-lineage");
	selectWindow("C1-b-lineage");
	setThreshold(1, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir4+"Lmx-b_"+title2+".tif");
	close("Lmx-b_"+title2+".tif");
	close("C1-b-lineage");
	
	/// Select a lineage nuclei	
	selectWindow("a-lineage");
	setTool("freehand");
	waitForUser("With FreeHand selection tool,\nselect a-lineage Lmx+ cells");
	
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
	
	run("Split Channels");
	close("C2-a-lineage");
	selectWindow("C1-a-lineage");
	setThreshold(1, 65535, "raw");
	run("Make Binary", "method=Default background=Default create");
	saveAs("tiff", dir5+"Lmx-a_"+title2+".tif");
	close("Lmx-a_"+title2+".tif");
	run("Close All");
	}
	} 
	else
	{ 
		//Create a binary file for Lmx
		print("Redo segmentation");
			

	}
	
}
	