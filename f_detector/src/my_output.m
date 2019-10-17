function []=my_output(header,detections, y_proc)
%Output Figures for each detection
%Plot with waveforms 
%--------------------------------------------------------------------------

%Create Directories
!mkdir -p WAVEFORMS

% define arrays
y_proc=cell2mat(y_proc);
% create plots : Processed waveforms
disp('Figures: Waveforms ..')
for k=1:length(detections)
fprintf('Figure: %03d out of %03d\n', k, length(detections))    
%----------------------------------------------------------
% Plot waveforms -- go through all available waveromfs
f1=figure('visible','off');
filename1=sprintf('waveforms_%03d.png',k);
for i=1:length(y_proc(1,:))
temp=y_proc(floor((detections(k,1)-5)./header(1).DELTA):round((detections(k,1)+15)./header(1).DELTA),i);
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
clear temp wvf %------------------------------------------------
end %end N of detections

end