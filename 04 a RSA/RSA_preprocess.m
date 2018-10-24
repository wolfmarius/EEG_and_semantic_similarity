
% load unique_words_wv.mat
% unique_words contains the list with all unique terms across all subjects (currently 225 for word2vec)
path_tmp = 'C:\data\marius\Code_and_Data\04 a RSA\';
load(strcat(path_tmp,  'unique_words_wv.mat'))
%load(strcat(path_tmp, 'unique_words_ngd.mat'))


% load the dissimilarity matrix. D is for NGD and D_w for word2ec.
% IMPORTANT: If you want to switch between D and D_wv you as well have
% switch line 20 and 21. --->  line20: [ac_D, ac_D_names] = cs_D(uniq_words, D);     % for D of NGD
% Moreover you as well have to switch the table that you load. line 39/40 load(strcat(path_table, selvps{n},'table_ngd.mat'))
%load(strcat(path_tmp, 'D_wv.mat'))
load(strcat(path_tmp, 'D_ngd_wv.mat'))   

%Switch as well the save function from ...wv.mat to ngd.mat if you change
%the input

% D is the similarity matrix for the normalized google distance. D_w is for word2vec
% [ac_D, ac_D_names] = ac_D(uniq_words, D_wv);       % for D_wv of w2v   
[ac_D, ac_D_names] = ac_D(uniq_words, D_ngd_wv);     % for D_ngd_wv of NGD


path_pre = 'C:\data\marius\02_RESULTS_sliding\';
path_out = 'C:\data\marius\04_RESULTS_rsa\';

path_table= 'C:\data\marius\02_RESULTS_sliding\tables\';

selvps   = {'vp01','vp02','vp03','vp04', 'vp05', 'vp06', 'vp07','vp08','vp09','vp10','vp11','vp12','vp13','vp14', 'vp15','vp16', 'vp17','vp18','vp19','vp20'};





for n=1:numel(selvps) 
 
    aaa = '_slidRSA_steps200ms_slide20ms_time0_1.mat'
    load(strcat(path_pre, selvps{n}, aaa))  
    load(strcat(path_table, selvps{n},'table_wv.mat'))
        
%% bring squiggle-codes and squiggle-names from EEG-Data in all_cat form 
clear all_cat
[all_cat(:,1),all_cat(:,2)]=find(ones(n_groups));


%all unique squiggle-codes
groups_code = unique(table1.trials_code);

%name of each squiggle
for i=1:size(groups_code, 1)
    for j=1:size(table1.trials_code,1)
    for k=2:size(table1.trials_code,1)
        if ((table1.trials_code(j) == groups_code(i)) && (~isequal(table1.trials_code(j), table1.trials_code(k))))
            
            groups_name(i) = table1.trials_name(j);    
        end           
    end
    end
end
groups_name = groups_name';       % contains the name of all consecutive categories in the table. e.g. 'Mensch', 'Auto', 'Katze, 'Mensch', 'Gesicht' ... 


%all_cat version of names and code
for a_ind=1:size(all_cat,1)
     a1=all_cat(a_ind,1);
     a2=all_cat(a_ind,2);
     
     ac_names{a_ind} = [groups_name{a1}, groups_name{a2},'-',groups_name{a2}, groups_name{a1} ];      
     ac_codes{a_ind} = [num2str(groups_code(a1)), num2str(groups_code(a2)),'-',num2str(groups_code(a2)), num2str(groups_code(a1))];         
end


% ac_names is necessary to extract the relevant squiggle-names from the
% dissimilarity matrix D so that you can get a dissimilarity matrix for
% each subject. ac_codes is not necessary in this analysis, but could be useful in a further analysis. 


%% Create D_vpN
 
% Bring D_vpN in corr_all-form  (this loop takes forever because of "contains")


Zeittest = 1
tic
for i=1:size(ac_names, 2)
    for j=1:size(ac_D_names, 2)    
        if contains(ac_names{i}, ac_D_names{j})  
           corr_all_D{i} =  ac_D(j);                   
        end
    end
end
toc
Zeittest = 2



%% Bring D_vpN from corr_all-form to corr_all_allsubs-form

corr_all_allsubs_D=zeros(size(t1,2), size(t2,2), size(corr_all_D, 2));
for i=1:size(corr_all_D, 2)      
    corr_all_allsubs_D(:,:,i) = repmat(corr_all_D{i}, size(t1,2)); 
end


%% Bring "corr_all" in "corr_all_allsubs"-form    (same process as in the permutation-script "permstat.m")
 
% [all_cat(:,1),all_cat(:,2)]=find(ones(n_groups));  % only excetute this line if you dont execute the first part of the script

 upper=ones(numel(t1));
        for ind=1:numel(t1)
           upper(:,ind)=[zeros(ind,1);ones(numel(t1)-ind,1)];
        end

        nan_ind=find(upper);
     
        [not_cat,ind_not] = setdiff(all_cat,cat_sparse,'rows');
        [in_cat,ind_in] = intersect(all_cat,cat_sparse,'rows');
   for i=1:size(all_cat,1)
        if isempty(find(ind_in==i))
           [x,ind_indata]=intersect(cat_sparse,[all_cat(i,2),all_cat(i,1)],'rows');
            corr_tmp=corr_all{ind_indata};
            %rotate and flip in order to switch upper diagonal with lower
            %diagonal
            corr_tmp=flip(rot90(corr_tmp));
            
        else
          [x,ind_indata]=intersect(cat_sparse,[all_cat(i,1),all_cat(i,2)],'rows');
            corr_tmp=corr_all{ind_indata};

        end

       % remove upper diagonal
       corr_tmp(nan_ind)=NaN;
       corr_all_allsubs(:,:,i)=corr_tmp;     
              
   end
   
   
   
%% save   
aaa = strcat(path_out,selvps{n},'_RSA_preprocess_steps',num2str(win*1000),'ms_slide',num2str(slide*1000),'ms_time',num2str(t_start),'_',num2str(t_end),'ngd225.mat')    
save(aaa,  'corr_all_allsubs', 'corr_all_allsubs_D','t1','t2', 'n_groups', 'groups_name', 'groups_code', 'ac_names', 'ac_codes', 'table1', 'win', 'slide', 't_start', 't_end', 'corr_all_D', 'aaa')  
      
%save(strcat(path_out,selvps{n},'_RSA_preprocess_steps',num2str(win*1000),'ms_slide',num2str(slide*1000),'ms_time',num2str(t_start),'_',num2str(t_end),'_ngd.mat'), 'corr_all_allsubs', 'corr_all_allsubs_D','t1','t2', 'n_groups', 'groups_name', 'groups_code', 'ac_names', 'ac_codes', 'table1', 'win', 'slide', 't_start', 't_end')   

%save(strcat(path_out,selvps{n},'_RSA_preprocess_steps',num2str(win*1000),'ms_slide',num2str(slide*1000),'ms_time',num2str(t_start),'_',num2str(t_end),'_ngd225.mat'), 'corr_all_allsubs', 'corr_all_allsubs_D','t1','t2', 'n_groups', 'groups_name', 'groups_code', 'ac_names', 'ac_codes', 'table1', 'win', 'slide', 't_start', 't_end')   


% "corr_all_allsubs" and "corr_all_allsubs_D" are the most important results.
% They are correlated and a permutation-test is performed on them in "RSA_perm_test.m"  
   
end
