//// This script requires the plugin 3Dsuite https://mcib3d.frama.io/3d-suite-imagej/#installation
//// Info on 3D suite: https://mcib3d.frama.io/3d-suite-imagej/uploads/neubias-3dsuite.pdf

setBatchMode(true); //This command operate the script without showing the images
dir1 = getDirectory("Image Folder");

index = lastIndexOf(dir1,"\\");
tempdir = substring(dir1,0,index);
dir2 = tempdir+"_LmxChecked/"; //directory of the Lmx binary files, changed named from _LmxCheckedBinary
dir3 = tempdir+"_MsxChecked/"; //directory of the Msx binary files, changed named from _MsxCheckedBinary

dir4 = tempdir+"_Stat/";
File.makeDirectory(dir4); //create thhe directory to save results
j = getFileList(dir1);
print("Number of images: "+j.length);
for (i=0; i<j.length; i++)
	{
	open(dir1 + j[i]);
	title1=getTitle();
	index1 = lastIndexOf(title1, ".");
	title2 =  substring(title1,0,index1);
	print(title2);
	
/////////// Measure Lmx intensity
	open(dir2 + "LmxBinary_"+title1);
	rename("Lmx1Binary");
	
	selectWindow(title1);
	run("Split Channels");
	close("C3-"+title2+".tif");
	
	selectWindow("C1-"+title2+".tif");
	rename("Lmx1");
	/// Substract background
	run("Duplicate...", "duplicate");
	rename("Lmx1 background");

	run("Gaussian Blur...", "sigma=35 stack");
	run("Multiply...", "value=1.25 stack");
	imageCalculator("Subtract create stack", "Lmx1","Lmx1 background");
	selectWindow("Result of Lmx1");
	rename("Lmx1Signal");
	
	/// Measure intensity using 3D Suite
	run("3D Intensity Measure", "objects=Lmx1Binary signal=Lmx1Signal");
		
	///////// Obtain total Lmx1 intensity
	Lmx1nuclei = nResults;
	print("Number of Lmx1+ nuclei: "+Lmx1nuclei);
	Lmx1intensity = 0;
		for (l = 0; l < nResults(); l++) 
		{
   		 v = getResult("IntensitySum", l);
   		 Lmx1intensity = Lmx1intensity+v;
		}
	print("Lmx1 Total intensity: "+Lmx1intensity);
	
	/// Add volume and centroid of each nuclei
	selectWindow("Lmx1Binary");
	run("3D Volume");
	run("3D Centroid");
	
	selectWindow("Results");
	saveAs("Results", dir4+"Lmx1Intensity_"+title2+".csv");
	run("Clear Results");
	close("Lmx1Signal");
	close("Lmx1Binary");
	close("Lmx1 background");
	
/////////// Measure Msx signal intensity
	open(dir3 + "MsxBinary_"+title1);
	rename("MsxBinary");
	selectWindow("C2-"+title2+".tif");
	rename("Msx");
	run("Duplicate...", "duplicate");
	rename("Msx background");
	run("Gaussian Blur...", "sigma=85 stack");
	imageCalculator("Subtract create stack", "Msx","Msx background");
	close("Msx");
	close("Msx background");
	selectWindow("Result of Msx");
	rename("MsxSignal");
	/// Measure intensity using 3D Suite
	run("3D Intensity Measure", "objects=MsxBinary signal=MsxSignal");
	
	///////// Obtain total Lmx1 intensity
	
	MsxIntensity = 0;
		for (k = 0; k < nResults(); k++) 
		{
   		 w = getResult("IntensitySum", k);
   		 MsxIntensity = MsxIntensity+w;
		}
	print("Msx Total intensity: "+MsxIntensity);
	/// Add volume and centroid of each "object"
	selectWindow("MsxBinary");
	run("3D Volume");
	run("3D Centroid");
	
	selectWindow("Results");
	saveAs("Results", dir4+"MsxIntensity_"+title2+".csv");
	run("Clear Results");
	run("Close All");
	
	}
	print("\nFinish");
	selectWindow("Log");
	saveAs("Text", dir4+"IntensitySummary_Log.csv");
	
setBatchMode(false);