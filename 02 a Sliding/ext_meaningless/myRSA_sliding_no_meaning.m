%% ADDING FIELDTRIP

%addpath F:\fieldtrip-20180108
addpath D:\Matlab_tools\fieldtrip-20170517
ft_defaults
%% LOADING DATA 
%path_pre= 'F:\Masterarbeit\data_marie_19.12\input_data\';
%path_out= 'F:\Masterarbeit\data_marie_19.12\scripts\1. Sliding\';
path_pre='C:\data\marius\readyforanalysis\';
path_behave = 'C:\data\marius\Code_and_Data\00 Preprocess_words\behave_blocks_w2v\';
path_out='C:\data\marius\02_RESULTS_sliding\';
mkdir(path_out)


%% VARIABLES FOR THE SLIDING WINDOW
t_start= 0;
t_end = 1;

win=0.2;   
slide=0.02; 

sr=50;

tois=t_start:1/sr:t_end;      %54  1 steps of 0.0050 from -0.3 to 2.4 if sr=200
t1=t_start:slide:(t_end-win); %lower bound for sliding window.  -0.3,  -0.29...2.35
t2=t1+win;                    %upper bound for sliding window.  -0.25, -0.24...2.40

ind_t1=1:slide/(1/sr):((numel(tois)-win/(1/sr))); %index for t1: step size depends on tois.  Roughly numel(tois)/numel(t1)
ind_t2=ind_t1+win/(1/sr);                         %index for t2: 


size_bins = win/(1/sr)+1;
n_bins=numel(t1);
n_channels = 64;

%%

%selvps   = {'vp03'};
%behave_b = {'vp03'};

selvps   = {'vp01', 'vp02','vp03','vp04', 'vp05', 'vp06', 'vp07','vp08','vp09','vp10', 'vp11','vp12','vp13','vp14', 'vp15','vp16', 'vp17', 'vp18','vp19', 'vp20'};
behave_b = {'vp01', 'vp02','vp03','vp04', 'vp05', 'vp06', 'vp07','vp08','vp09','vp10', 'vp11','vp12','vp13','vp14', 'vp15','vp16', 'vp17', 'vp18','vp19', 'vp20'};


for n=1:numel(selvps)

load(strcat(path_pre, selvps{n}  , '_cleaned'));
load(strcat(path_behave, behave_b{n}, '_k'));    
    
%% Extraction Meaningful
vpn = selvps{n};
[table1_meaningless] = ext_meaningless(behave_blocks, data, vpn);              

'

% Extract and sort trials 
    %Sort by code    
    table1_meaningless = sortrows(table1_meaningless,'trials_code','ascend');
    code_groups = findgroups(table1_meaningless.trials_code);
    code_blocks = table1_meaningless.trials_block;

    %Sort by name 
    % table1 = sortrows(table1,'trials_name','ascend')
    % name_groups = findgroups(table1.trials_name);

%Extract trials by name or by code:
trials_index = table1_meaningless{:,2};


%code_groups = code_groups(1:20);
%trials_index_reduced = trials_index(1:20);

%% PREPROCESSING


   cfg=[];
cfg.resamplefs=sr;
cfg.detrend='no';
data=ft_resampledata(cfg,data);


cfg=[];
cfg.latency=[t_start t_end];
cfg.trials = trials_index;
data=ft_selectdata(cfg, data);
      
cfg=[];
cfg.keeptrials='yes';
%cfg.keeptrials='no';
data=ft_timelockanalysis(cfg, data);

n_trials=size(data.trialinfo,1);

%z-trans across trials 
mean_trials=repmat(nanmean(data.trial,1),n_trials,1,1);      
std_trials=repmat(nanstd(data.trial,1),n_trials,1,1);  
data.trial=(data.trial-mean_trials)./std_trials;

%figure()
%plot(data.time,squeeze(nanmean(data.trial))');  %this plots the ERP of a
%single person

 data_vec=zeros(n_trials, n_bins, numel(data.label)*size_bins);  %creates the vector in which you save everything.
    for bin=1:n_bins
       % vectorize sel_bins: data_vec(trials, nbins, features)
       data_vec_tmp=data.trial(:,:,ind_t1(bin):ind_t2(bin));
       data_vec(:,bin,:)=reshape(data_vec_tmp,n_trials,[]);         % trials*n_bins* [st_ep(bin_size)*channels]=64 channels in each bin. 
                                     
    end

% cat_sparse
   clear corr_mat_cat cat1 cat2 cat_sparse corr_all auto_trial
   n_groups = max(code_groups);
   corr_mat_cat=zeros(n_groups);
        for ind=1:n_groups % this loop creates a square matrix with 1 on the diagonal and in the upper triangle
           corr_mat_cat(:,ind)=[ones(ind,1);zeros(n_groups-ind,1)];
        end
         [cat1,cat2]=find(corr_mat_cat); %returns row and column index of non-zero elements
          cat_sparse=[cat1,cat2];
          
          
for c_ind=1:size(cat_sparse,1) %for loop executed with parallel computing.    
               [corr_all{c_ind}, auto_trial{c_ind}]=vectorized_sliding_after_no_meaning(data, data_vec, c_ind, n_bins, cat_sparse, code_groups, code_blocks);  %kann ich durch eine korrelation ersetzen.               
end





save(strcat(path_out,selvps{n},'_slidRSA_steps',num2str(win*1000),'ms_slide',num2str(slide*1000),'ms_time',num2str(t_start),'_',num2str(t_end),'_nomeaning.mat'), 'tois','corr_all','cat_sparse','auto_trial','t1','t2','n_groups','code_groups', 'win', 'slide', 't_start', 't_end')

end
