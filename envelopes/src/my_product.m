function     trace=my_product(sub_env)
%Calclualte product of envelopes for
% each subarray
%--------------------------------------------------------------------------
array=cell2mat(sub_env);
array=array';

trace=prod(array); %product

trace=trace'; %transpose

end