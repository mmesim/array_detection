function   [sub_env]=my_envelope(stack)
%function to calculate normalized envelopes 
%for each subarray


%Preallocate memory
sub_env=cell(1,length(stack));

parfor i=1:length(stack)
temp=stack{1,i};
temp_trans=hilbert(temp);
temp_env=abs(temp_trans);

sub_env{1,i}=temp_env./max(temp_env);
    
end
end