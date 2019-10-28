function [Fcrit,Fscaled,c]=my_fstatistics(F,lcorner,hcorner,T,header)
%Function to compute Fstatistic
%N1=2BT [B: bandwidth, T: window length [sec]]
%N2=2BT(N-1) [N:Number of stations]
%and search for normalized factor [c]
%see Arrowsmith et al., 2009 [BSSA] & Arrowsmith et al., 2008 [GJI]
%--------------------------------------------------------------------------
delta=header(1).DELTA;  %sampling interval
N=length(header);        %N of stations

%Get bandwidth
if lcorner==hcorner
    B=(1./(2*delta))-lcorner;  % F-Nyqist - lower cut off
else
    B=hcorner-lcorner;   % high - low
end
%Degrees of Freedom
N1=2*B*T; N2=N1*(N-1);

%Theoretical Distribution
x=0:0.01:10;
y=fpdf(x,N1,N2);
[~,ind]=max(y);
xmax=x(ind);

%Critical value at 99% F(0.99)
Fcrit=finv(0.99,N1,N2);

%Search for C coefficient
c=1;
fmax=0;
while abs(fmax-xmax)>0.2 
c=c-0.0001;
%---------- If c<0 then change line 32 -------------------------------
if c<0 
     msg='Error with c calculations! Check my_fstatistics.m line 32.';
     error(msg)
end
%---------------------------------------------------------------------
% Count normalized F and find maximum value
[n,edges]=histcounts(c.*F);
[~,ind]=max(n);
fmax=(edges(ind+1)+edges(ind))./2;
end


%F normalized by c
Fscaled=c.*F;


end
