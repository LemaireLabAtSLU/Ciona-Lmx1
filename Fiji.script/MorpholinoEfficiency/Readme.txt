This directory contains the Fiji script to quantify signal intensity from the fluorescent reporter Lmx1 with Morpholino site>H2B::neonGreen and Msx>mCherry.
First, z-stacks of neurula embryos (2048x2048 pixels) were saved from czi files to tiff with the macro "1.1.MacroSaveTiff_halfPixelResolution_fromCzi.ijm", which reduces the resolution by half (1024x1024pixels).
The H2B::NG signal and the mCherry signals were then binarized using the macro "2.Binary.ijm". 
The binary was check, any noise from outside the endogenous dorsal domain of Msx and Lmx1 expression was removed with the macro "3.Check.Binary.ijm".
The intensity of the neonGreen and mCherry signal was measured using the script "4.IntensityMeasurement3D.ijm".
The total neongreen intensity was then divided by the mCherry total intensity per embryos to account for variation of the injected volume.
