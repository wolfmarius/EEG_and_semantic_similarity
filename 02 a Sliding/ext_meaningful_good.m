function [trials_index, table1, trials_index_nomeaning] = ext_meaningful_good(behave_blocks, data, vpn)


NGD=behave_blocks(:,9);

for i=1:size(behave_blocks,1)
       if isequal(behave_blocks{i,9}, 'n') 
           index_n(i) = i;
       end
end

ind_named_bb = setdiff(1:size(NGD), index_n);

b = behave_blocks(:,1);


squiggle_code = cell2mat(b(ind_named_bb)); % code of all  named squiggles
squiggle_name =NGD(ind_named_bb);          % name of all  named squiggles             

%% Index of all named squiggles in DATA
trials_code_all  = data.trialinfo(:,8); 

trials_index= find(ismember(trials_code_all, squiggle_code)); %index of trials with a named squiggle in data.trialinfo
trials_code = trials_code_all(ismember(trials_code_all, squiggle_code));  %code of all named trials in DATA


% All named trials
%trials_all_named = data.trial(trials_index)';         %all named trials

%% Name of all named trials
b1 = cell2mat(behave_blocks(:,1));
b9 = behave_blocks(:,9); 

k = 1;
n_of_trials = length(trials_code);
tmp_name={zeros(n_of_trials,1)}';
for j=1:n_of_trials
    for i=1:120 
            if trials_code(j)==b1(i)
             tmp_name(k) = b9(i);
              k = k+1;
            end
    end
end
trials_name = tmp_name';  %name of the named trials in DATA

%Blocks
block_tmp = data.trialinfo(:,7);
trials_block= block_tmp(trials_index(:));


%name_code = horzcat(num2cell(trials_index), trials_name, num2cell(trials_code));
vpn = repmat(vpn, length(trials_name),1);
table1 = table(vpn, trials_index, trials_name, trials_code, trials_block);   

%% Index of meaningless trials
trials_index_nomeaning = setdiff(1:max(trials_index), trials_index)';




end

