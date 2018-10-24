%this script executes ext_meaningful_good.m to extract the meaningful
%trials and saves them for each
%subject into a table. 
%The tables are loaded in myRSA_sliding_table_already.m


path_pre = 'C:\data\marius\Code_and_Data\00 Preprocess_words\behave_blocks_w2v\'  ;   % for w2v
%path_pre = 'C:\data\marius\Code_and_Data\00 Preprocess_words\behave_blocks_NGD\'  ;   % for NGD

path_data='C:\data\marius\readyforanalysis\';

path_out='C:\data\marius\02_RESULTS_sliding\';

selvps   = {'vp01', 'vp02','vp03','vp04', 'vp05', 'vp06', 'vp07','vp08','vp09','vp10', 'vp11','vp12','vp13','vp14', 'vp15','vp16', 'vp17', 'vp18','vp19', 'vp20'};
behave_b = {'vp01', 'vp02','vp03','vp04', 'vp05', 'vp06', 'vp07','vp08','vp09','vp10', 'vp11','vp12','vp13','vp14', 'vp15','vp16', 'vp17', 'vp18','vp19', 'vp20'};
%selvps   = {'vp18'};
%behave_b = {'vp18'};

for n=1:numel(selvps)

load(strcat(path_data, selvps{n}  , '_cleaned'));
load(strcat(path_pre, behave_b{n}, '_k'));    
    
%% Extraction Meaningful
vpn = selvps{n};
[trials_index, table1, trials_index_nomeaning] = ext_meaningful_good(behave_blocks, data, vpn);              

% Extract Meaningless
    %rials_index_meaningless = setdiff(1:size(data.trial,2), trials_index)'

% Extract and sort trials 
    %Sort by code    
    table1 = sortrows(table1,'trials_code','ascend');
    code_groups = findgroups(table1.trials_code);
    code_blocks = table1.trials_block;

    %Sort by name 
    % table1 = sortrows(table1,'trials_name','ascend')
    % name_groups = findgroups(table1.trials_name);

%Extract trials by name or by code:
trials_index = table1{:,2};


save(strcat(path_out,selvps{n},'table_wv.mat'), 'table1')
%save(strcat(path_out,selvps{n},'table_ngd.mat'), 'table1')

end