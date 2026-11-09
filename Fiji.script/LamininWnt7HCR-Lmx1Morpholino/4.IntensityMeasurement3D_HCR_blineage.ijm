//// This script requires the plugin 3Dsuite https://mcib3d.frama.io/3d-suite-imagej/#installation
//// Info on 3D suite: https://mcib3d.frama.io/3d-suite-imagej/uploads/neubias-3dsuite.pdf

//setBatchMode(true); //This command operate the script without showing the images


dir1 = getDirectory("Image Folder");

index = lastIndexOf(dir1,"\\");
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_Lmx_blin/"; //directory of the b10.73-76 cell binary files, 
dir3 = tempdir+"_Lmx_alin/"; //directory of the a11.93-96 cell binary files, 

dir4 = tempdir+"_blin-LamininStat/";
File.makeDirectory(dir4); //create the directory to save results for b nuclei - laminin signal
dir5 = tempdir+"_blin-Wnt7Stat/";
File.makeDirectory(dir5); //create the directory to save results for b nuclei - Wnt7 signal
dir6 = tempdir+"_blin-LamininSignal/";
File.makeDirectory(dir6); //create the directory to save object Laminin signal images in b nuclei 
dir7 = tempdir+"_blin-Wnt7Signal/";
File.makeDirectory(dir7); //create the directory to save object Wnt7 signal images in b nuclei 

j = getFileList(dir1);
print("Number of images: "+j.length);


run("3D Manager Options", "volume integrated_density mean_grey_value std_dev_grey_value mode_grey_value minimum_grey_value maximum_grey_value distance_between_centers=10 distance_max_contact=1.80 drawing=Contour");// run the manager 3D and add image
run("3D Manager");
Ext.Manager3D_LiveRoi();

for (i=0; i<j.length; i++)
	{

		
	open(dir1 + j[i]);
	title1=getTitle();
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
		
	/// Find background for Laminin and Wnt7
	Stack.setActiveChannels("10110");
	Stack.setPosition(1,1,1); 
	run("Specify...", "width=100 height=100 x=256 y=256 slice=1");
 	waitForUser("select background area"); //place area in trunk 
	run("Clear Results");
	run("Measure");
	Wnt7Bg = getResult("Mean", 0);
	print("Wnt7 Background: "+Wnt7Bg);
	d=getSliceNumber();
	Stack.setPosition(4,d,1); 
	run("Measure");
	LamininBg = getResult("Mean", 1);
	print("Laminin Background: "+LamininBg);
	run("Clear Results");
	
	selectWindow(title1);
	run("Split Channels");
	close("C3-"+title2+".tif");
	close("C2-"+title2+".tif");
	close("C5-"+title2+".tif");
	selectWindow("C1-"+title2+".tif");
	rename("Wnt7");
	selectWindow("C4-"+title2+".tif");
	rename("Laminin");	
	
	
	/// Substract background for laminin
	selectWindow("Laminin");
	run("Subtract...", "value="+LamininBg+" stack");
	rename("LamininSignal");
	
	setMinAndMax(0, 2000);
	run("8-bit");
	//close("Laminin background");
	//close("Laminin");
	
	/// Substract background for Wnt7
	selectWindow("Wnt7");
	run("Subtract...", "value="+Wnt7Bg+" stack");
	rename("Wnt7Signal");
	setMinAndMax(0, 2000);
	run("8-bit");
	close("Wnt7 background");
	close("Wnt7");
	
/////////// Measure  nascent intensity in b lineage
	print("b nuclei measurements");
	open(dir2 + "Lmx-b_"+title1);
	rename("Lmx1bBinary");
	
	///// Select only nuclei
	imageCalculator("AND create stack", "LamininSignal","Lmx1bBinary");
	close("LamininSignal");
	selectWindow("Result of LamininSignal");
	rename("LamininSignal");
	
	imageCalculator("AND create stack", "Wnt7Signal","Lmx1bBinary");
	close("Wnt7Signal");
	selectWindow("Result of Wnt7Signal");
	rename("Wnt7Signal");
	
/// Check that there are nuclei

selectWindow("Lmx1bBinary");
run("Z Project...", "projection=[Max Intensity]");
selectWindow("MAX_Lmx1bBinary");
setThreshold(128, 255);
run("Create Selection");
run("Measure");
LmxMax = getResult("Max", 0); /// if the image is empty, it is selecting all the image anyway and the maximum pixel intesity will be 0
print("Lmx binary area max intensity: "+ LmxMax);

if (LmxMax == 0)
{
print("No Lmx-positive nuclei");
run("Close All");
run("Clear Results");
}
else 
{

close("MAX_Lmx1bBinary");
run("Clear Results");
	
	///Transfer 3D ROI
selectWindow("Lmx1bBinary");
run("3D Objects Counter", "threshold=128 slice=17 min.=100 max.=146800640 objects statistics summary");	

selectWindow("Statistics for Lmx1bBinary");
Table.rename("Statistics for Lmx1bBinary", "Results");
bnuclei = nResults; //// Find number of object
print(bnuclei);

selectWindow("Results");
saveAs("Results", dir6+title2+"Lmx1-bNuclei.csv");
run("Close");


selectImage("Objects map of Lmx1bBinary");
Ext.Manager3D_AddImage;
Ext.Manager3D_Count(nb);
print("Number of object:"+nb);

selectImage("Objects map of Lmx1bBinary");
run("glasbey_inverted");
saveAs("tiff", dir6+title2+"_b-nuclei.tif");
selectWindow(title2+"_b-nuclei.tif");
run("Close");

//////Measure nascent transcription for individual a-nucleus 

///////Create individual mask nuclei
///Need to remove signal outside of nuclei of selection by creating mask of individual nuclei

selectImage("Lmx1bBinary");
    
getDimensions(width, height, channels, slices, frames);

for (p = 0; p<bnuclei; p++) 
{
	print("Nucleus #"+p+1);
	
	newImage("Nucleus-"+p+1, "8-bit black", width, height, slices);

	// if list is not visible please refresh list by using Deselect
	Ext.Manager3D_Select(0);
    Ext.Manager3D_DeselectAll();
    
    Ext.Manager3D_Select(p);
   	Ext.Manager3D_FillStack(255, 255, 255);
    
	setThreshold(128, 255);
	run("Convert to Mask", " ");

}




//////Measure Laminin intensity
    print("Laminin intensity measurement");
   
for (s = 0; s<bnuclei; s++) 
{
print("Nucleus #"+s+1);

	selectImage("LamininSignal");
    run("Duplicate...", "duplicate");
	rename("LamininSignalNucleus-"+s+1);

	imageCalculator("AND create stack", "LamininSignalNucleus-"+s+1, "Nucleus-"+s+1);
	close("LamininSignalNucleus-"+s+1);

	selectWindow("Result of LamininSignalNucleus-"+s+1);
	rename("LamininSignalNucleus-"+s+1);
	
	//crop image, if not most thresholding are set to 0
	selectWindow("LamininSignalNucleus-"+s+1);
	Ext.Manager3D_Select(0);
    Ext.Manager3D_DeselectAll();
    Ext.Manager3D_Select(s);
    run("Crop");
	
	run("Duplicate...", "duplicate");
	rename("Threshold");

	///Threshold for spot
	resetThreshold();
	selectWindow("Threshold");
	run("Auto Local Threshold", "method=Bernsen radius=30 parameter_1=0 parameter_2=0 white stack");
	
	////select appropriate volume
	run("3D Objects Counter", "threshold=255 slice=12 min.=25 max.=15000 objects statistics summary");
	selectWindow("Statistics for Threshold");
	Table.rename("Statistics for Threshold", "Results");
	lamininNascentSpotb = nResults;
	saveAs("Results", dir4+"LamininSpotVolume-"+title2+"-Nucleus-"+s+1+".csv");
	run("Close");
	print("Number of Laminin Spot: "+lamininNascentSpotb);
	
	close("Threshold");
	selectWindow("Objects map of Threshold");
	rename("Objects");
	
	///if result == 0, no nascent transcription
if(lamininNascentSpotb == 0) 
{
	print("No Nascent transcription");
}	
//// else use object as mask for 3D intensity measurements
	else 
{	
	run("3D Intensity Measure", "objects=Objects signal=LamininSignalNucleus-"+s+1);
	
	for (q = 0; q < lamininNascentSpotb; q++) 
	{
	Nascent = getResult("IntensitySum", q);
	print("Laminin nascent spot-"+q+1+" of Nucleus #"+s+1+": "+Nascent);
	}
	
	saveAs("Results", dir4+"LamininSpotIntensity-"+title2+"-Nucleus-"+s+1+".csv");
	run("Close");
}	
	selectWindow("Objects");
	run("glasbey_inverted");
	saveAs("tiff", dir6+title2+"-NascentLamininSignalNucleus-"+s+1+"-map.tif");
	selectWindow(title2+"-NascentLamininSignalNucleus-"+s+1+"-map.tif");
	run("Close");
	selectWindow("LamininSignalNucleus-"+s+1);
	saveAs("tiff", dir6+title2+"-NascentLamininSignalNuclei-"+s+1+".tif");
	selectWindow(title2+"-NascentLamininSignalNuclei-"+s+1+".tif");
	run("Close");
}

	 //////Measure Wnt7 intensity
    print("Wnt7 intensity measurement");

for (r = 0; r<bnuclei; r++) 
{
	print("Nucleus #"+r+1);
	selectImage("Wnt7Signal");
	
		//select only one nucleus
	run("Duplicate...", "duplicate");
	rename("Wnt7SignalNucleus-"+r+1);
	selectImage("Wnt7SignalNucleus-"+r+1);
	imageCalculator("AND create stack", "Wnt7SignalNucleus-"+r+1, "Nucleus-"+r+1);
	close("Wnt7SignalNucleus-"+r+1);
	selectWindow("Result of Wnt7SignalNucleus-"+r+1);
	rename("Wnt7SignalNucleus-"+r+1);
	
	//crop image, if not most thresholding are set to 0
	selectWindow("Wnt7SignalNucleus-"+r+1);
	Ext.Manager3D_Select(0);
    Ext.Manager3D_DeselectAll();
    Ext.Manager3D_Select(r);
    run("Crop");

	run("Duplicate...", "duplicate");
	rename("Threshold");

	///Threshold for spot
	resetThreshold();
	selectWindow("Threshold");
	run("Auto Local Threshold", "method=Bernsen radius=30 parameter_1=0 parameter_2=0 white stack");
	
	////select appropriate volume
	run("3D Objects Counter", "threshold=255 slice=12 min.=25 max.=15000 objects statistics summary");
	selectWindow("Statistics for Threshold");
	Table.rename("Statistics for Threshold", "Results");
	Wnt7NascentSpotb = nResults;
	saveAs("Results", dir5+"Wnt7SpotVolume-"+title2+"-Nucleus-"+r+1+".csv");
	run("Close");
	print("Number of Wnt7 Spot: "+Wnt7NascentSpotb);
	
	close("Threshold");
	selectWindow("Objects map of Threshold");
	rename("Objects");
	
	///if result == 0, no nascent transcription
if(Wnt7NascentSpotb == 0) 
{
	print("No Nascent transcription");
}	
//// else use object as mask for 3D intensity measurements
	else 
{	
	run("3D Intensity Measure", "objects=Objects signal=Wnt7SignalNucleus-"+r+1);
	
	for (k = 0; k < Wnt7NascentSpotb; k++) 
	{
	Nascent = getResult("IntensitySum", k);
	print("Wnt7 nascent spot-"+k+1+" of Nucleus #"+r+1+": "+Nascent);
	}


	saveAs("Results", dir5+"Wnt7SpotIntensity-"+title2+"Nuclei-"+r+1+".csv");
	run("Close");
}
	selectWindow("Objects");
	run("glasbey_inverted");
	saveAs("tiff", dir7+title2+"-NascentWnt7SignalNuclei-"+r+1+"-map.tif");
	selectWindow(title2+"-NascentWnt7SignalNuclei-"+r+1+"-map.tif");
	run("Close");
	selectWindow("Wnt7SignalNucleus-"+r+1);
	saveAs("tiff", dir7+title2+"-NascentWnt7SignalNucleus-"+r+1+".tif");
	selectWindow(title2+"-NascentWnt7SignalNucleus-"+r+1+".tif");
	run("Close");


	
}
}
	Ext.Manager3D_Reset(); ///Remove nuclei from ROI manager

	run("Clear Results");
	run("Close All");
	
	}
	
	
	print("\nFinish");
	selectWindow("Log");
	saveAs("Text", dir6+"IntensitySummary_blineage_Log.csv");
	
//setBatchMode(false);