function   [Sx,Sy,F,baz,S,N] = my_beam(header, Sxmax, Sxinc, Sxmin, Symax, Symin, Syinc,rlon,rlat,y_proc,win)
%Create a beam for each grid point using Nth root stacking
%--------------------------------------------------------------------------
% Create grid with all possible pairs
[Gx,Gy] = meshgrid(Sxmin:Sxinc:Sxmax, Symin:Syinc:Symax);
pairs = [Gx(:),Gy(:)];
% -------------------------------
%Preallocate memory
stlo=zeros(length(header),1);
stla=zeros(length(header),1);
%extract longitude and latitude for each station
for i=1:length(header)
stlo(i,1)=header(i).STLO;
stla(i,1)=header(i).STLA;
end
%---------------------------------------------------------------
y=cell2mat(y_proc); % from cell to array;   
%define window to cut waveform
mwindow=round(win./header(1).DELTA);
%-----------------------------------------------------------------
n=length(y)-mwindow;
N=fix(1:mwindow/2:n);

parfor ii=1:length(N)
%Calculate Fmax for each chunk
[Sx(ii,1),Sy(ii,1),F(ii,1),baz(ii,1),S(ii,1)]=my_shifts(stlo,stla,rlon,rlat,pairs,y(N(ii):N(ii)+mwindow,:),header(1).DELTA);

end

end