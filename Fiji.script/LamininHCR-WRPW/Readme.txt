This directory contains the Fiji script to quantify Laminin HCR signal intensity.
First, z-stacks of mid-tailbud embryos were "cleaned" by removing any additional embryos from the stack. Dorsal neural tube nuclei expressing Lmx1 reporter were then segmented using “nuclei_3d_binary.ijm”. 
The segmentation was then manually checked and corrected. 
The total intensity of the Laminin signal was measured using the script "Intensity measurement_based on binary.ijm". 
It was then divided by the number of segmented nuclei counted with "Intensity measurement_based on binary.ijm" to obtain the average Laminin intensity per nucleus for each embryo.
