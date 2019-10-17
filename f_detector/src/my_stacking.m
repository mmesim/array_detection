function beam=my_stacking(ynew,waveform,N)
%calculate beam using F detector (Blandford et al., 1974)
%-----------------------------------------------------------
%fix waveform length
ynew=ynew(1:length(waveform(1,:)),:)';

%N is the number of stations
a=(N-1)./N;


%Sum shifted waveforms
b=(sum(waveform)).^2;
%residuals between shifted and sum of input
c=1/N*(sum(ynew));

d=sum((waveform-c).^2);

%F beam 
beam=a.*(b./d);

end