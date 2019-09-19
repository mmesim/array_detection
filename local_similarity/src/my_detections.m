function  [detections]=my_detections(network_trace,det_win,thres,delta,time_thres) 
%Declare detections when amlitude>MD+(thres*MAD)

%define window
my_win=round(det_win./delta);

%Preallocate memory
N=length(network_trace(:,1))-my_win;
ind=cell(N,1);

%------- Start Detections -------------------------------------------------
parfor i=1:N

my_window=network_trace(i:my_win+i,1);
MAD=mad(my_window,1); % Median absolute deviation
MD=median(my_window); % median 

threshold=MD+(thres*MAD); % detection threshold (amplitude)

index=find(my_window(:,1)>=threshold); % Apply detection threshold  

ind{i,1}=[index+i ones(length(index),1)*MAD ones(length(index),1)*MD my_window(index)];

end
%--------------------------------------------------------------------------
% Group detections -- time threshold defined by user
% See time_thres variable in parameter.m

%find unique values and calculate N time MAD
table=cell2mat(ind);
table=[table (table(:,4)-table(:,3))./table(:,2)];
[unq,ic,ia]=unique(table(:,1));
%------------------------------------------------------------
mad_times=table(ic,5);
%convert to sec and group detections
det=(unq*delta);
inter=diff(det); %interevent times
ind3=find(inter(:,1)>time_thres); %apply threshold
%--------------------------------------------------------------------------
det2=[det(1,1); det(ind3(2:end),1)];  %in seconds

%maximum MAD times
j=1;
for k=1:length(ind3)
mad_times2(k,1)=max(mad_times(j:ind3(k)));    
j=ind3(k);    
%mad_times2=[mad_times(1,1); mad_times(ind3(2:end),1)];

end
%final array
detections=[det2 mad_times2];
end
