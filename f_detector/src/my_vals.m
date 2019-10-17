function [baz,S]=my_vals(Sx,Sy)
%Calculate back azimuth and slowness vector
%--------------------------------------------------------------------------

%Slowness vector
S=sqrt(Sx.^2+Sy.^2);

%Backazimuth
if Sx<0 && Sy<0
baz=90-57.295*(pi/2);

elseif Sx<0 && Sy>0
baz=90-57.295*(3.0*pi/2);

else 

baz=90-57.296*atan2(Sy,Sx);
end

if baz<0
baz=baz+360;
end

end
