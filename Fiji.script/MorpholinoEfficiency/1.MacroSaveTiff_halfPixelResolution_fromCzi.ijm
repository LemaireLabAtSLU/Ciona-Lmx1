dir1 = getDirectory("Choose Source Directory");
j=nImages;
print("number of images: "+j);
for (i=0; i<j; i++)
{
title1=getTitle();
	index2 = lastIndexOf(title1,".");
	title2 = substring(title1,0,index2);

getDimensions(w, h, channels, slices, frames); 
run("Size...", "width="+w/2+" height="+h/2+" depth="+slices+" constrain average interpolation=Bilinear");
	saveAs("tiff", dir1+title2+".tif");
close();
}
