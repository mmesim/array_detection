function []=my_catalog(detections,header)
%Output catalog using header information and datetime function
%--------------------------------------------------------------------------
time=datetime(header(1).KZTIME,'Format','HH:mm:ss.SSS'); %header time
t=time+seconds(detections(:,1)); %time in seconds since the beginning 
%---------------------------------------------------------
%output catalo
fout=fopen('detections.txt','w');
for i=1:length(t)
fprintf(fout,'%s %s \n', header(1).KZDATE ,t(i));
end
fclose(fout);

end