function  [my_stack]=my_stacking(y, ind, header)
%Devide array into subarrays and perform linear stacking
%without time shifting
%Stations are sorted by their number as defined in parameter.m 
%ind variable
%This will work for nodals that usually have namess 1:N.
%--------------------------------------------------------------------------

%Sort Stations by station number
for i=1:length(header)
stations(i,1)=str2double(header(i).KSTNM);
end
[~,index]=sort(stations);
%use ind to sort the waveforms by station number
y=y(1,index);

%preallocate
my_stack=cell(1,length(ind(1,:)));
temp=cell2mat(y);

for i=1:length(ind(1,:))
%Define subarrays    
ind2=ind{:,i};   
wv=temp(:,ind2);
%Normalize each waveforms
for k=1:length(wv(1,:))
wv2(:,k)=wv(:,k)./max(wv(:,k));
end

wv2=wv2';
%Stack
stack=(sum(wv2))';
%Normalize stacked waveforms
my_stack{1,i}=stack./max(stack);

clear wv2    
end


end