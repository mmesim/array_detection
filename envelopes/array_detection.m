% Package to detect events using envelopes of stacked sub-arrays          % 
% The code is based on Meng and Ben-Zion 2018 (GJI)                       %
%                                                                         %
% --------------- M.Mesimeri 07/2019 --------------------------------------
clear;clc;close all; tic %start timer

%% 00.Setup
parameters %load parameter file
parpool('local',workers); %Start parallel pool
mydir=pwd; pdir=sprintf('%s/src/',pwd); % get working directory path
addpath(genpath(pdir)); %add all *.m scripts to path
%--------------------------------------------------------------------------
%% 01. Load data (Sac files)
display('Loading files..')
[y,header]=my_loadfiles(mydata);

%% 02. Preprocess data 
display('Preprocessing...')
y_proc=my_preprocessing(y,header(1).DELTA,type,lcorner,hcorner);

%% 03. Stack sub-arrays
display('Stacking...')
stack=my_stacking(y_proc,ind,header);

%% 04. Envelopes
display('Envelopes...')
sub_env=my_envelope(stack);

%% 05. Evelopes' product
display('Product...')
trace=my_product(sub_env);

%% 06. Detect events using STA_LTA
display('Detecting events...')
[detections,ratio]=my_sta_lta(trace,sta,lta,thres,time_thres,header(1).DELTA);

%% 07. Output #1 : Catalog 
display('Output #1 : Catalog')
my_catalog(detections,header);

%% 08. Output #2 : Figures
display('Output #2 : Figures')
saved_wvf=my_output(header,detections, y_proc,ratio,thres);

%% 09. Shutdown parallel pool
delete(gcp)
fprintf('Elapsed times %6.2f minutes... \n',toc/60) %stop timer



