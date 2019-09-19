function [saved_wvf]=my_output(header,detections, y_proc,ratio,thres)
%Output Figures for each detection
%(1) Plot with waveforms 
%(2) Plot with sta/lta ratio for each detection
%(3) Sta/lta
%Output saved_wvf can be used later for beamforming
%--------------------------------------------------------------------------

%Create Directories
mkdir WAVEFORMS
mkdir STA_LTA

% define arrays
y_proc=cell2mat(y_proc);
% create plots : Processed waveforms, sta/lta 
disp('Figures: Waveforms, and sta/lta ..')
for k=1:length(detections)
fprintf('Figure: %03d out of %03d\n', k, length(detections))    
%----------------------------------------------------------
% Plot waveforms -- go through all available waveromfs
f1=figure('visible','off');
filename1=sprintf('waveforms_%03d.png',k);
for i=1:length(y_proc(1,:))
temp=y_proc(floor((detections(k,1)-5)./header(1).DELTA):round((detections(k,1)+15)./header(1).DELTA),i);
saved_wvf{k,1}{i,1}=temp;
wvf=temp./max(temp);
plot((1:length(wvf))*header(1).DELTA,wvf+i,'k');
hold on
end
hold off
xlim([0 length(wvf)*header(1).DELTA])
ylim([1 length(header)])
xlabel('Time [sec]','FontSize',16)
ylabel('Station ID','FontSize',16)
saveas(f1,filename1,'png');
movefile(filename1,'WAVEFORMS/');
close all
%-------------------------------------------------
%%STA/LTA ratio
f2=figure('visible','off');
filename2=sprintf('sta_lta%03d.png',k);
temp2=ratio(floor((detections(k,1)-5)./header(1).DELTA):round((detections(k,1)+15)./header(1).DELTA));
plot((1:length(wvf))*header(1).DELTA, temp2,'k');
hold on
xlim([0 length(wvf)*header(1).DELTA])
ylim([min(temp2) max(temp2)])
xlabel('Time [sec]','FontSize',16)
ylabel('STA/LTA','FontSize',16)
%--------------------------------------------------
saveas(f2,filename2,'png');
movefile(filename2,'STA_LTA/');
close all
clear loc temp wvf temp2
%------------------------------------------------
end %end N of detections
% STA/LTA -- one figure
f3=figure('visible','off');
filename3=sprintf('sta_lta_%d_%d.png',header(1).NZYEAR,header(1).NZJDAY);
plot((1:length(ratio))*header(1).DELTA, ratio);
hold on
N=thres*median(ratio);
X=1:length(ratio)*header(1).DELTA; Y=N*ones(length(X),1);
plot(X,Y,'m:','LineWidth',2)
title(sprintf('Date:%d-%d Median STA/LTA: %5.3f',header(1).NZYEAR,header(1).NZJDAY,median(ratio)))
xlim([0 length(ratio)*header(1).DELTA])
ylim([min(ratio) max(ratio)])
xlabel('Time [sec]','FontSize',16)
ylabel('STA/LTA','FontSize',16)
set(gca,'FontSize',16)
%--------------------------------------------------
saveas(f3,filename3,'png');
close all

end