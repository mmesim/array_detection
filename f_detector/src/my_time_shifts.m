function [tau] = my_time_shifts(stlo,stla,rlon,rlat,kx,ky)
%01. Calculate time shifts

%Azimuth
az=azimuth(rlat,rlon,stla,stlo);

%Distance from reference point
%dist=sqrt((x-xo).^2+(y-yo).^2);
dist=distance(stla,stlo,rlat,rlon);

% Calculate rx and ry   
rx=dist*(cosd(90-az));
ry=dist*(sind(90-az));

%Get time shift in seconds
tau=((kx*rx)+(ky*ry));


end