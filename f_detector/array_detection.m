% Package to detect events using F-detector                               % 
% The code is based on Blandford (1974,Geophysics)                        %
% see also Arrowsmith et al. (2008 GJI;2009 BSSA)                         %                                                               %
% --------------- M.Mesimeri 10/2019 --------------------------------------
clear;clc;close all; tic %start timer

%% 00.Setup
parameters %load parameter file
parpool('local',workers); %Start parallel pool
mydir=pwd; pdir=sprintf('%s/src/',pwd); % get working directory path
addpath(genpath(pdir)); %add all *.m scripts to path
%--------------------------------------------------------------------------
%% 01. Load data (Sac files)
disp('Loading files..')
[y,header]=my_loadfiles(mydata);

%% 02. Preprocess data 
disp('Preprocessing...')
y_proc=my_preprocessing(y,header(1).DELTA,type,lcorner,hcorner);

%% 03. F detector
disp('Beam...')
[Sx,Sy,F,baz,S,time] = my_beam(header, Sxmax, Sxinc, Sxmin, Symax, Symin, Syinc,rlon,rlat,y_proc,win);

%% 04. F statistics
disp('F statistics...')
[Fcrit,Fscaled,c]=my_fstatistics(F,lcorner,hcorner,win,header);
%% 05. Detections
disp('Detections...')
[detections,table]=my_detections(Fscaled,Fcrit,time,time_thres,Sx,Sy,baz,S,header);
%% 06. Output #1 : Catalog
disp('Output #1 : Catalog...')
my_catalog(detections,header);

%% 07. Output #2 - Figures
disp('Output #2 : Figures')
my_output(header,detections, y_proc);

%% 08. Shutdown parallel pool
delete(gcp)
fprintf('Elapsed time %6.2f minutes... \n',toc/60) %stop timer
save results.mat


