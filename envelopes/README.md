
## Parameter File
**mydata**: string, full path to waveforms. The directory should contain all the waveforms for one day in SAC format. If you want to loop through different days you can edit the code or run it multiple times. 

**workers**: integer, define the number of cores for parallel processing.                 

**ind**: array to define different subarrays based on the station number. For nodals station name is a number.  

**type**: string, choose between 'low', 'high', 'bandpass' for filtering, for high or low pass set the hcorner equal to lcorner and vice versa

**lcorner**: lower corner frequency

**hcorner**: higher corner frequency

**sta**: Define short window in seconds

**lta**: Define long term window in seconds

**thres**:   define the detection threshold (thres * median(sta/lta) )

**time_thres**: define time threshold to group detections in sec


*Package written in Matlab for event detection using Large N array*

The algorithm is based on: 
* Meng, H., Ben-Zion, Y., 2018. Detection of small earthquakes with dense array data: Example from the San Jacinto fault zone, southern California. Geophys. J. Int. 212, 442–457. https://doi.org/10.1093/gji/ggx404 *

## Execute
Edit the parameter file. Double check that the path is correct. Type "array_detection" in the command window. 

parameter.m, array_detection.m & src should be in the same directory. Waveforms can be anywhere.



