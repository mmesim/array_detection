%Parameter file for array_detection code                                   %
%                                                                          % 
%---------- M.Mesimeri 07/2019 --------------------------------------------

%path to waveforms
mydata='/uufs/chpc.utah.edu/common/home/koper-group1/keith/YL_microseisms_2018/Data_nodal/2018.06.28/04'; 
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
sta=0.1;            % Short term window [sec]
lta=2;              % Long term window [sec]
thres=10;           % Detection threshold (thres * median sta/lta)
time_thres=1;       % Time threshold to group detections[sec]
%--------------------------------------------------------------------------
%Output file : Define window around the detection
out_win=20;          %Window around the detection in sec
