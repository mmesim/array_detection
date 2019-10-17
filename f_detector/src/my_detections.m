function  [detections,table]=my_detections(Fscaled,Fcrit,time,time_thres,Sx,Sy,baz,S,header)
%function to group detections and associate phases
%returns a vector with detections in seconds 
%--------------------------------------------------------------------------
delta=header(1).DELTA;
%find detections above F(0.99)
index=find(Fscaled(:,1)>=Fcrit);

%Correct time and set it equal to the center of each sliding window
for i=1:length(time)-1
ntime(i,1)=(time(1,i+1)+time(1,i))./2; 
end
ntime=[ntime ; (time(1,end)+75)./2];

det=(ntime(index).*delta);

Sx1=Sx(index); Sy1=Sy(index); baz1=baz(index); S1=S(index);

%--------------------------------------------------------------------------
%Group detections
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
%return table with Sx Sy baz [backazimuth] S [Slowness vector)
table=[Sx1(ind) Sy1(ind) baz1(ind) S1(ind)];

end