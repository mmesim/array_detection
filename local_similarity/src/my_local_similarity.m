function loc_sim=my_local_similarity(ind,y_proc,cc_win,delta)
%calculate local similarity between nearest neighbors

%Preallocate memory
loc_sim=cell(length(ind(1,:)),1);

%Running windows
my_win=round(cc_win./delta);

for i=1:length(ind(1,:))
% grab waveforms for each group
fprintf('Group_%03d \n',i)
x=y_proc(ind(:,i));
%------------------------------------------------    
%Correlations with the MASTER station (always x{1,1})
loc_sim{i,1}=my_correlation(my_win,x);
end

end