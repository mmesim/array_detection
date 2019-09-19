function [y_proc]=my_preprocessing(y,delta,type,lcorner,hcorner)
%Preprocessing
%(1) remove mean
%(2) remove trend
%(3) filter
%(4) fix waveform length
%--------------------------------------------------------------------------
%Preallocate memory
y_proc1=cell(1,length(y));
y_fix=cell(1,length(y));


parfor i=1:length(y)
temp=y{1,i};
%01. remove mean
ym=temp-mean(temp);
%02. remove trend
yr=my_detrend(ym,1);
%03. filter
yf=my_filter(yr,type,delta,lcorner,hcorner);
y_proc1{1,i}=yf;
end

%fix waveform length
val=min(cellfun('length',y_proc1));

for k=1:length(y)
y_fix{1,k}=y_proc1{1,k}(1:val,:);   
 
end

y_proc=y_fix;
end