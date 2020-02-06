# Multi Frame Super Resolution Tool

A MFSR Tool for Matlab inspired by the "Robust And Fast Super Resolution" tool by Oded Hanson [Link](http://www1.idc.ac.il/toky/videoProc-07/projects/SuperRes/srproject.html)

![MFSR_Tool](Images/MFSR_Tool.png)

## Info

Tool to compute an image of higher resolution from a video of low resolution images.  
Choose from multiple Image-Registration methods and Super-Resolution algorithms.

Supported video input formats:

- .avi
- .mov
- .mp4
- .m4v

## Source Code

The source code is located in `MFSR_App` folder.

## Standalone Apps

Compiled standalone apps for Windows and OS X are located in `MFSR_Tool` folder.  
If Matlab (R2019b) isn't installed, use `MFSR_Installer_web` to download the Matlab runtime.

## Video Samples

Some video samples can be found in `Video_Samples` to test the different algorithms  
Most of them are downloaded from:

- MDSP Super-Resolution And Demosaicing Datasets [Link](https://users.soe.ucsc.edu/~milanfar/software/sr-datasets.html)

## Sources

All Papers of the used algorithms can be found in the folder `Papers` or here:

- Video super-resolution by adaptive kernel regression [Link](https://www.mathworks.com/matlabcentral/fileexchange/60766-video-super-resolution-by-adaptive-kernel-regression)
- Fast and Robust Multiframe Super Resolution [Link](http://people.duke.edu/~sf59/SRfinal.pdf)
- Pyramidal Implementation of the Lucas Kanade Feature Tracker Description of the algorithm [Link](http://robots.stanford.edu/cs223b04/algo_tracking.pdf)
