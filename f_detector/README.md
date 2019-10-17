## Parameter File
**mydata**: string, full path to waveforms. The directory should contain all the waveforms for one day in SAC format. If you want to loop through different days you can edit the code or run it multiple times.

**workers**: integer, define the number of cores for parallel processing.

**type**: string, choose between 'low', 'high', 'bandpass' for filtering, for high or low pass set the hcorner equal to lcorner and vice versa

**lcorner**: lower corner frequency

**hcorner**: higher corner frequency

**rlon**: Longitude of the reference point (in Decimal degrees)

**rlat**: Latitude of the reference point (in Decimal degrees)

**win**: sliding window for detection

**Sxmin**: Minimun value for Sx

**Sxmax**: Maximum value for Sx

**Sxinc**: Increment for Sx

**Symin**: Minimun value for Sy

**Symax**: Maximum value for Sy

**Syinc**: Increment for Sy

**time_thres**: define time threshold to group detections in sec

## Execute
Edit the parameter file. Double check that the path is correct. Type "array_detection" in the command window.

parameter.m, array_detection.m & src should be in the same directory. Waveforms can be anywhere. 

The code by default uses a 50% overlapping window. 
