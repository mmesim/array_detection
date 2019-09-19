function network_trace=my_stacking(loc_sim)
%convert cell to matrix and sum
%there is no need to normalize, values are already [0 1]

%Play with loc_sim array
loc_sim=loc_sim';
loc_sim=cell2mat(loc_sim);
loc_sim=loc_sim';

%Stack
stack=(sum(loc_sim))';
stack=stack./max(stack);

%--------------------------------------------
%remove long period trend
%1st order polynomial

network_trace=my_detrend(stack,1);

end 