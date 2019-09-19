function [ind]=my_clustering(neighb,header)
%Calcualte distance between all stations and 
%choose the N nearest neighbors
%It's easy to replace this function with distance function
%--------------------------------------------------------------------------

%Preallocate memory
stlo=zeros(length(header),1);
stla=zeros(length(header),1);
ind=zeros(neighb+1,length(stla));

%extract longitude and latitude for each station
for i=1:length(header)
stlo(i,1)=header(i).STLO;
stla(i,1)=header(i).STLA;
end

%find nearest neighbors
parfor k=1:length(stla)
%Master Station    
xo=(stlo(k,1).*(cos(stla(k,1))*pi/180))*111.11;
yo=stla(k,1)*111.11;   

%All stations
x=(stlo.*(cos(stla)*pi/180))*111.11;
y=stla*111.11; 

%Distance from master station
dis=sqrt((x-xo).^2+(y-yo).^2);
%Sort distances and keep the first three
%First is 0 (Master station)
%Second and third are the two nearest neighbors
[~,index]=sort(dis);
%save in to a matrix
ind(:,k)=index(1:neighb+1); 
end
end
