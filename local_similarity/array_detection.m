% Package to detect events using local similarity for Large arrays        % 
% The code is based on the Li et al., 2018 (Scientific Reports)           %
%                                                                         %
% --------------- M.Mesimeri 07/2019 --------------------------------------
clear;clc;close all; tic %start timer

%% 00.Setup
parameters %load parameter file
parpool('local',workers); %Start parallel pool
mydir=pwd; pdir=sprintf('%s/src/',pwd); % get working directory path
addpath(genpath(pdir)); %add all *.m scripts to path
%-----------------------------htop---------------------------------------------
%% 01. Load data (Sac files)
display('Loading files..')
[y,header]=my_loadfiles(mydata);

%% 02. Preprocess data 
display('Preprocessing...')
y_proc=my_preprocessing(y,header(1).DELTA,type,lcorner,hcorner);

%% 03. Nearest neighbor and master station groups
display('Clustering...')
ind=my_clustering(neighb,header);
s
%% 04. Perform cross correlations and define local similarity for each group
display('Local similarity...');
loc_sim=my_local_similarity(ind,y_proc,cc_win,header(1).DELTA);

%% 05. Stack all the local similarity outputs
display('Stacking...')
network_trace=my_stacking(loc_sim);

%% 06. Detect events  above a threshold level
display('Detecting events...')
detections=my_detections(network_trace,det_win,thres,header(1).DELTA,time_thres);

%% 07. Output #1 -- Catalog
display('Output #1 : Catalog')
my_catalog(detections,header)

%% 08. Output #2 -- Figures
display('Output #2 : Figures')
[saved_loc,saved_wvf]=my_output(header,detections, y_proc,loc_sim,network_trace);
%% 09. Shutdown parallel pool
delete(gcp)
fprintf('Elapsed time %6.2f minutes... \n',toc/60) %stop timer

clearvars -except detections saved_loc saved_wvf header
save detections.mat
