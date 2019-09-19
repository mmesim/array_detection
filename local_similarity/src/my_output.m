function [saved_loc,saved_wvf]=my_output(header,detections, y_proc,loc_sim,network_trace)
%Output Figures for each detection
%(1) Image of local similarity 
%(2) Plot with waveforms
%(3) Plot with network trace
%Output asved_loc and save_wvf can be used later for beamforming
%--------------------------------------------------------------------------

%% Create Directories
mkdir WAVEFORMS
mkdir LOCAL_SIMILARITY
mkdir NETWORK_TRACE

%% Sort Stations by station number
for i=1:length(header)
stations(i,1)=str2double(header(i).KSTNM);
end
[sorted_sta,ind]=sort(stations);
%use ind to sort the waveforms by station number
y_proc=y_proc(1,ind);
loc_sim=loc_sim(ind,1);
%% define arrays
loc_sim=cell2mat(loc_sim'); 
y_proc=cell2mat(y_proc);
%% create plots : Processed waveforms, local similarity, and network trace
display('Figures: Waveforms, Local similarity, and Network Trace..')
for k=1:length(detections)
fprintf('Figure: %03d out of %03d\n', k, length(detections))    
loc=loc_sim(floor((detections(k,1)-5)./header(1).DELTA):round((detections(k,1)+15)./header(1).DELTA),:);
loc=loc';
saved_loc{k,1}=loc;
%%-------- image of local similarity
%--------------------------------------------------------------------------
f1=figure('visible','off');
filename1=sprintf('local_similarity_%03d.png',k);
image(loc,'CDataMapping','scaled');
colormap(flipud(gray));
ticks=0:1250:5000;
labels=ticks*header(1).DELTA;
set(gca, 'XTick', ticks, 'XTickLabel', labels);
xlabel('Time [sec]','FontSize',16)
ylabel('Station ID','FontSize',16)
colorbar
set(gca,'YDir','normal')
saveas(f1,filename1,'png');
movefile(filename1,'LOCAL_SIMILARITY/');
close all
%----------------------------------------------------------
% Plot waveforms -- go through all available waveromfs
f2=figure('visible','off');
filename2=sprintf('waveforms_%03d.png',k);
for i=1:length(y_proc(1,:))
temp=y_proc(floor((detections(k,1)-5)./header(1).DELTA):round((detections(k,1)+15)./header(1).DELTA),i);
saved_wvf{k,1}{i,1}=temp;
wvf=temp./max(temp);
plot((1:length(wvf))*header(1).DELTA,wvf+i,'k');
hold on
end
hold off
xlim([0 length(wvf)*header(1).DELTA])
ylim([1 length(stations)])
xlabel('Time [sec]','FontSize',16)
ylabel('Station ID','FontSize',16)
saveas(f2,filename2,'png');
movefile(filename2,'WAVEFORMS/');
close all
%-------------------------------------------------
%% Network trace
f3=figure('visible','off');
filename3=sprintf('network_trace%03d.png',k);
temp2=network_trace(floor((detections(k,1)-5)./header(1).DELTA):round((detections(k,1)+15)./header(1).DELTA));
plot((1:length(wvf))*header(1).DELTA, temp2,'k');
xlim([0 length(wvf)*header(1).DELTA])
ylim([min(temp2) max(temp2)])
xlabel('Time [sec]','FontSize',16)
ylabel('Local coherence','FontSize',16)
%--------------------------------------------------
saveas(f3,filename3,'png');
movefile(filename3,'NETWORK_TRACE/');
%close all
clear loc temp wvf temp2
%------------------------------------------------
end %end N of detections

end