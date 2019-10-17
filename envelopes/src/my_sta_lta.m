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
j=1;
for ii=1:length(det)-1
if det(ii+1)-det(ii)<time_thres
temp(j)=ii+1;
j=j+1;
end
end
ind=setdiff(1:length(det)-1,temp);

%Final vector with detections in seconds
detections=det(ind);

end
