function  [cc]=my_correlation(my_win,x)
%correlations bewtween two stations

%length of parfor
N=length(x{1,1})-my_win;

%Preallocate memory
cc2=zeros(N,1); 

for k=1:length(x)-1

parfor i=1:N
% %GPU use     
% %----Send to GPU
% X=gpuArray(x{1,1}(i:my_win+i,1));
% Y=gpuArray(x{1,k+1}(i:my_win+i,1));
% %-Do xcorr
% [cc1,~]=xcorr(X,Y,'coeff'); 
% %--- Bring back to cpu
% r=gather(cc1);
% cc2(i,k)=max(abs(r));

%CPU use
[cc1,~]=xcorr(x{1,1}(i:my_win+i,1),x{1,k+1}(i:my_win+i,1),'coeff');    
cc2(i,k)=max(abs(cc1));
end

end

%Mean CC for the traces -- local similarity 
cc=(mean(cc2'))';
end
