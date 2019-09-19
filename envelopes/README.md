
## Parameter File
**mydata**: string, full path to waveforms. The directory should contain all the waveforms for one day in SAC format. If you want to loop through different days you can edit the code or run it multiple times. 

**workers**: integer, define the number of cores for parallel processing.                 

**type**: string, choose between 'low', 'high', 'bandpass' for filtering, for high or low pass set the hcorner equal to lcorner and vice versa

**lcorner**: lower corner frequency

**hcorner**: higher corner frequency

**neighb**: define the number of nearest neighbors

**cc_win**: define the window for cross correlation in sec 

**thres**:   define the detection threshold (median+(thres*MAD))

**det_win**: define the detection window in sec
  
**time_thres**: define time threshold to group detections in sec


*Package written in Matlab for event detection using Large N array*

The algorithm is based on: 
*Li, Z., Peng, Z., Hollis, D., Zhu, L., McClellan, J., 2018. High-resolution seismic event detection using local similarity for Large-N arrays. Sci. Rep. 8, 1–10. https://doi.org/10.1038/s41598-018-19728-w

## Execute
Edit the parameter file. Double check that the path is correct. Type "array_detection" in the command window. 

parameter.m, array_detection.m & src should be at the same directory. Waveforms can be anywhere.



