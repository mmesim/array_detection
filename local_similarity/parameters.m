%Parameter file for array_detetion code                                   %
%                                                                         % 
%---------- M.Mesimeri 06/2019 --------------------------------------------

%path to waveforms
mydata='/uufs/chpc.utah.edu/common/home/koper-group1/mesimeri/Array_Detection/LOCAL_SIMILARITY/DATA1'; 
%--------------------------------------------------------------------------
% Parallel settings
workers=32;                 %Set number of cores to work on a local machine
%-------------- Filtering parameters --------------------------------------
type='high';                %'low', 'high', 'bandpass'
lcorner=1;                  % lower corner frequency
hcorner=1;                  % higher corner frequency
%----------------- Clustering parameters ----------------------------------
neighb=4;                   % Number of nearest neighbors
%----------------- Correlation Parameters ---------------------------------
cc_win=2;                   % Window for cross correlation [sec]
%----------------- Detection Parameters -----------------------------------
thres=8;                    % Detection threshold (median+(thres*MAD))
det_win=60;                 % Detection window [sec]
time_thres=15;              % Time threshold to group detections [sec]
%--------------------------------------------------------------------------
