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
[unq,ic,~]=unique(table(:,1));
%------------------------------------------------------------
mad_times=table(ic,5);
%convert to sec and group detections
det=(unq*delta);

j=1;
for ii=1:length(det)-1
if det(ii+1)-det(ii)<time_thres
temp(j)=ii+1;
j=j+1;
end
end
ind2=setdiff(1:length(det)-1,temp);

%Final vector with detections in seconds
det=det(ind2);


%maximum MAD times
j=1;
for k=1:length(ind2)
mad_times2(k,1)=max(mad_times(j:ind2(k)));
j=ind2(k);
%mad_times2=[mad_times(1,1); mad_times(ind3(2:end),1)];

end
%final array
detections=[det mad_times2];
end

