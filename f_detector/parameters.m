%Parameter file for array_detetion code                                   %
%                                                                         % 
%---------- M.Mesimeri 10/2019 --------------------------------------------

%path to waveforms
mydata='DATA1'; 
%--------------------------------------------------------------------------
% Parallel settings
workers=32;                 %Set number of cores to work on a local machine
%-------------- Filtering parameters --------------------------------------
type='high';                %'low', 'high', 'bandpass'
lcorner=1;                  % lower corner frequency
hcorner=1;                  % higher corner frequency
%--------------------------------------------------------------------------
% Reference Point (Center of the Array)
rlon =-112.887; rlat=38.500;
% Define sliding window in seconds
win=3; 
% Search bounds for Sx and Sy (sec/deg)
Sxmin=-20; Sxmax=20; Sxinc=1;
Symin=-20; Symax=20; Syinc=1;
%--------------------------------------------------------------------------
% Detection Parameters 
time_thres=15;              % Time threshold to group detections [sec]
%--------------------------------------------------------------------------
