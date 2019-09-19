%Parameter file for array_detection code                                   %
%                                                                          % 
%---------- M.Mesimeri 07/2019 --------------------------------------------

%path to waveforms
mydata='/uufs/chpc.utah.edu/common/home/koper-group1/mesimeri/Array_Detection/LOCAL_SIMILARITY/DATA3'; 
%--------------------------------------------------------------------------
% Parallel settings
workers=32;                  %Set number of cores to work on a local machine

%--------------------------------------------------------------------------
%set groups
ind={[1:50]' [56:65 151]' [51:55 66:80]' [81:95]' [96:110]' [111:126]' ...
    [127:140]' [141:145]' [146:150]'}; 
%-------------- Filtering parameters --------------------------------------
type='high';    %'low', 'high', 'bandpass'

lcorner=1;           % lower corner frequency
hcorner=1;           % higher corner frequency
%----------------- Detection Parameters -----------------------------------
sta=1;               % Short term window [sec]
lta=10;              % Long term window [sec]
thres=7;             % Detection threshold (thres * median sta/lta)
time_thres=15;       % Time threshold to group detections[sec]
%--------------------------------------------------------------------------
