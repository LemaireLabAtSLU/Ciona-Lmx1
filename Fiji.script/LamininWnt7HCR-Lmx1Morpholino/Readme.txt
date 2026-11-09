This directory contains the Fiji script to quantify Laminin and Wnt7 HCR signal intensity on initial tailbud embryos injected with Lmx1 morpholino resistant>H2B::NG and either Lmx1 morpholino or control morpholino.
First, z-stacks of initial-tailbud embryos (2048x2048 pixels) were saved from czi files to tiff with the macro "1.MacroSaveTiff.ijm".
The H2B::NG signal corresponding to dorsal neural tube nuclei was then binarized using the macro "2.Binary.ijm". 
The nuclei were then assigned as derived from the a- or b-lineage with the macro "3.Check.Binary.ijm".
The intensity of the Laminin and Wnt7  nuclear spots was measured using the script "4.IntensityMeasurement3D_HCR_blineage.ijm" for the b-derivered dorsal neural tube cells. 
The two brightest spots per nuclei were then sum (total in one nucleus). It was then summed for the whole embryo and divided by the number of nuclei. The signal was then normalized to the mean of the control for each btach (in percent).

