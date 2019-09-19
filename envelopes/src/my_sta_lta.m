function [detections,ratio]=my_sta_lta(trace,sta,lta,thres,time_thres,delta)
%Run sta/lta detector for the envelopes
%Declare detections when sta/lta>thres*median(sta/lta)
%thres is defined in parameter.m 
%--------------------------------------------------------------------------
%Set windows to samples
sta=fix(sta/delta);
lta=fix(lta/delta);

%Pre allocate memory
N=length(trace);
ratio=zeros(N,1);
%------------------------------------------------------
%Run sta/lta -- Attention: trace is already an evelope
   parfor i=lta+1:N
       mlta=mean(trace(i-lta:i));
       msta=mean(trace(i-sta:i));
       ratio(i)=msta/mlta;
   end
%Find Potential triggers
triggs=find(ratio(:,1)>=thres*median(ratio));

%Group triggers 
det=(triggs*delta);
inter=diff(det); %interevent times
ind=find(inter(:,1)>time_thres); %apply threshold

%Final vector with detections in seconds
detections=det(ind);
   
end