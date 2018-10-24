function [table1_nomeaning] = ext_meaningless(behave_blocks, data, vpn)

%this function is executed in  myRSA_sliding_no_meaning.m 
%It extracts the indices of all meaningless trials and saves them into a
%table for each subject


NGD=behave_blocks(:,9);

for i=1:size(behave_blocks,1)
       if isequal(behave_blocks{i,9}, 'n') 
           index_n(i) = i;
       end
end

ind_named_bb = setdiff(1:size(NGD), index_n); %hier evt. 'n' hinzufügen
index_no_meaning_bb = setdiff(1:size(NGD), ind_named_bb);
b = behave_blocks(:,1);


squiggle_code_nomeaning = cell2mat(b(index_no_meaning_bb)); % code of all  meaningless squiggles


%% Index of all named squiggles in DATA
trials_code_all  = data.trialinfo(:,8); 

trials_index= find(ismember(trials_code_all, squiggle_code_nomeaning)); %index of trials with no meaning squiggle in data.trialinfo
trials_code = trials_code_all(ismember(trials_code_all, squiggle_code_nomeaning));  %code of all meaningless trials


%Blocks
block_tmp = data.trialinfo(:,7);
trials_block= block_tmp(trials_index(:));


%name_code = horzcat(num2cell(trials_index), trials_name, num2cell(trials_code));
vpn = repmat(vpn, length(trials_index),1);


table1_nomeaning = table(vpn, trials_index, vpn, trials_code, trials_block);   


%% Index of meaningless trials





end

