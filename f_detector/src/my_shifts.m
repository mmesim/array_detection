function [Sx,Sy,F,baz,S]=my_shifts(stlo,stla,rlon,rlat,pairs,y,delta)
%Calculate time shifts for each grid point
%Shift waveforms using fft
%Return beam after stacking using Nth root stacking
%--------------------------------------------------------------------------

for ii=1:length(pairs)
    
for i=1:length(stlo)
    
%1. Get time shifts
tau=my_time_shifts(stlo(i),stla(i),rlon,rlat,pairs(ii,1),pairs(ii,2));

%2. Shift waveforms 
waveform(i,:)=my_fft(tau,y(:,i),delta);

end
%3.Stacking 
beam(ii,:)=my_stacking(y,waveform,length(stlo));

end
%4.Maximum Amplitude
[Sx,Sy,F] = my_ampls(beam,pairs);
%5.Baz and Slowness
[baz,S]=my_vals(Sx,Sy);
end
