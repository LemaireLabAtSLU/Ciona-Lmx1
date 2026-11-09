dir1 = getDirectory("Choose Source Directory");
j=nImages;
for (i=0; i<j; i++)
{
title1=getTitle();
	index2 = lastIndexOf(title1,".");
	title2 = substring(title1,0,index2);
saveAs("tiff", dir1+title2+".tif");
close();
}
