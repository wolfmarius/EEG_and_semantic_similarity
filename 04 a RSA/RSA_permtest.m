% Permutation test for RSA

addpath D:\Matlab_tools\fieldtrip-20170517
ft_defaults

%%

nrand=1000;
alpha=0.05;
alpha_c=0.05;

path_pre = 'C:\data\marius\04_RESULTS_rsa\';
path_out = 'C:\data\marius\04_RESULTS_rsa\';

selvps   = {'vp01', 'vp02','vp03','vp04', 'vp05', 'vp06', 'vp07','vp08','vp09','vp10', 'vp11','vp12','vp13','vp14', 'vp15','vp16', 'vp17', 'vp18','vp19', 'vp20'};


for n=1:numel(selvps)
aaa =     '_RSA_preprocess_steps200ms_slide20ms_time0_1_wv.mat'
load(strcat(path_pre, selvps{n},aaa))

corr_all_allsubs_tmp   = permute(corr_all_allsubs,[3,1,2]);   %corr_all_allsubs is the temporal generalization matrix for one subject
corr_all_allsubs_D_tmp = permute(corr_all_allsubs_D,[3,1,2]); %corr_all_allsubs_D is the semantic dissimialrity matrix for one subject


corr_all_allsubs_perm{n}     =  corr_all_allsubs_tmp;
corr_all_allsubs_D_perm{n}  =   corr_all_allsubs_D_tmp;

end

corr_all_allsubs    = corr_all_allsubs_perm;
corr_all_allsubs_D  = corr_all_allsubs_D_perm;


%% RSA with permutation

tic
for r=1:nrand   %does the randomisation
    
for n=1:numel(selvps)

[d1,d2,d3] = size(corr_all_allsubs_perm{n});
idx = randperm(d1);
corr_all_allsubs_D_perm{n} = corr_all_allsubs_D_perm{n}(idx,:,:);   %shuffles the first dimension of the corr_all_allsubs_D_perm


for i=1:size(t1,2)
    for j=1:size(t1,2)          
        [RHO_perm{n}(i,j), PVAL_perm{n}(i,j)] = corr(corr_all_allsubs_perm{n}(:, i, j), corr_all_allsubs_D_perm{n}(:, i, j), 'type','Spearman', 'rows', 'complete');     
        tval_perm(n,i,j) = RHO_perm{n}(i,j) / sqrt((1-RHO_perm{n}(i,j)^2) / (size(corr_all_allsubs_D_perm{n}, 1) - 2)); % Transform RHO into t-value to get a t-map (see Bortz p.223) 
                       
        %[TAU{n}(i,j), PVAL_TAU{n}(i,j)] = corr(corr_all_allsubs{n}(:, i, j), corr_all_allsubs_D{n}(:, i, j), 'type','Kendall', 'rows', 'complete');
        %z{n}(i,j)      = norminv(PVAL_TAU{n}(i,j)) * sign(TAU{n}(i,j));
        % %the test statistic is approximatley normal distributed. Thats why
        % %we can transform them into z-values (Bortz p. 224)      
    end
end
end


tval_perm = tval_perm*(-1);
[h,p,c, stat] = ttest(tval_perm, 0,'Alpha', alpha_c);  



   
      mask_neg=squeeze(h.*(stat.tstat<0));
   rep=find(isnan(mask_neg));
   mask_neg(rep)=0;
   [L_neg,num_neg] = bwlabel(mask_neg); 
    for neg=1:num_neg
        m=find(L_neg==neg);
        negt(neg)=sum(stat.tstat(m));
    end
    if num_neg==0
        negt=0;
    else
    end
    
     mask_pos=squeeze(h.*(stat.tstat>0));
     mask_pos(rep)=0;
    [L_pos,num_pos] = bwlabel(mask_pos);   
    for pos=1:num_pos
        m=find(L_pos==pos);
        post(pos)=sum(stat.tstat(m));
    end
    
     if num_pos==0
        post=0;
    else
    end
        
    pos_tsum{r}=sort(post,'descend');
    neg_tsum{r}=sort(negt,'ascend');
    
    clear post negt  
      
end
 toc

%% find clusters

min_pos=min(cellfun(@numel,pos_tsum));
min_neg=min(cellfun(@numel,neg_tsum));

for x=1:nrand
   pos_tsum{x}= pos_tsum{x}(1:min_pos);
   neg_tsum{x}= neg_tsum{x}(1:min_neg);
end

pos_dist=reshape([pos_tsum{:}],min_pos,nrand);
pos_dist=sort(pos_dist,2,'descend');

neg_dist=reshape([neg_tsum{:}],min_neg,nrand);
neg_dist=sort(neg_dist,2,'ascend');


%% RSA without permutation: Correlate "corr_all_allsubs" and "corr_all_allsubs_D"

for n=1:numel(selvps)


for i=1:size(t1,2)
    for j=1:size(t1,2)
        [RHO{n}(i,j), PVAL{n}(i,j)] = corr(corr_all_allsubs{n}(:, i, j), corr_all_allsubs_D{n}(:, i, j), 'type','Spearman', 'rows', 'complete');     
        tval(n,i,j) = RHO{n}(i,j) / sqrt((1-RHO{n}(i,j)^2) / (size(corr_all_allsubs_D{n}, 1) - 2)); % Transform RHO into t-value to get a t-map (see Bortz p.223)         
                
        %[TAU{n}(i,j), PVAL_TAU{n}(i,j)] = corr(corr_all_allsubs{n}(:, i, j), corr_all_allsubs_D{n}(:, i, j), 'type','Kendall', 'rows', 'complete');
        %z{n}(i,j)      = norminv(PVAL_TAU{n}(i,j)) * sign(TAU{n}(i,j));
        % %the test statistic is approximatley normal distributed. Thats why
        % %we can transform them into z-values (Bortz p. 224)      
                   
    end
end
end

%%
tval = tval*(-1);
[h,p,c, stat] = ttest(tval, 0,'Alpha', alpha_c);


    
    
    
        data_mask_neg=squeeze(h.*(stat.tstat<0));
        data_mask_neg(rep)=0;     
        [data_L_neg,data_num_neg] = bwlabel(data_mask_neg);
        if data_num_neg==0
           data_num_neg=1;
        end           
    for neg=1:data_num_neg
        m=find(data_L_neg==neg);
        data_negt(neg)=sum(stat.tstat(m));
    end
        [data_negt,ind_negt]=sort(data_negt,'ascend');
    
        
        data_mask_pos=squeeze(h.*(stat.tstat>0));
        data_mask_pos(rep)=0;     
        [data_L_pos,data_num_pos] = bwlabel(data_mask_pos);    
        if data_num_pos==0
           data_num_pos=1;
        end    
    for pos=1:data_num_pos
        m=find(data_L_pos==pos);
        data_post(pos)=sum(stat.tstat(m));
    end
       [data_post,ind_post]=sort(data_post,'descend');
 
       
       
     mask_pos=squeeze(h.*(stat.tstat>0));
     mask_pos(rep)=0;
    [L_pos,num_pos] = bwlabel(mask_pos);  
    for pos=1:num_pos
        m=find(L_pos==pos);
        post(pos)=sum(stat.tstat(m));
    end
    
     if num_pos==0
        post=0;
    else
    end
                   
        
        % test stepwise for significance        
       i=0;
       sig=1;
       sig_poscluster=0;
       mask_pos=zeros(size(data_mask_pos));   
while sig==1 i<min_pos
       i=i+1;     
       sig= data_post(i)>=pos_dist(1,round(alpha*nrand));
       p_pos(i)=nearest(pos_dist,data_post(i))/nrand;
       
       if sig==1
       sig_poscluster(i)=i;
       mask_pos=mask_pos+(data_L_pos==ind_post(i));
       end
end
              
        
        i=0;
        sig=1;
        sig_negcluster=0;
        mask_neg=zeros(size(data_mask_neg));
while sig==1 & i<min_neg
       i=i+1;     
       sig= data_negt(i)<=neg_dist(1,round(alpha*nrand));
       p_neg(i)=nearest(neg_dist,data_negt(i))/nrand;

       if sig==1
       sig_negcluster(i)=i;
       mask_neg=mask_neg+(data_L_neg==ind_negt(i));
       end
end   
        
       
%         i=0;
%         sig=1;
%         sig_poscluster=0;
%         mask_pos=zeros(size(data_mask_pos));        
% while sig==1 i<min_pos
%        i=i+1;     
%        sig= data_post(i)>=pos_dist(i,round(alpha*nrand*(1/i)));
%        
%        if sig==1
%        sig_poscluster(i)=i;
%        mask_pos=mask_pos+(data_L_pos==ind_post(i));
%        end
% end
%        
%         
%         
%         i=0;
%         sig=1;
%         sig_negcluster=0;
%         mask_neg=zeros(size(data_mask_neg));
% while sig==1 & i<min_neg
%        i=i+1;     
%        sig= data_negt(i)<=neg_dist(i,round(alpha*nrand*(1/i)));
%       if sig==1
%       sig_negcluster(i)=i;
%       mask_neg=mask_neg+(data_L_neg==ind_negt(i));
%       end
% end   

%% Plot
        mask_alpha=mask_pos+mask_neg;
        ind_z=find(mask_alpha==0);
        mask_alpha(ind_z)=0.5;
        
        %%
        figure(1)
   H= imagesc((t1+t2)*0.5,(t1+t2)*0.5,squeeze(stat.tstat),[-3 3]); 
       set(gca,'YDir','normal')
       colorbar
       set(H,'AlphaData',mask_alpha) 
       title('NGD,  sampling rate=100HZ,  frequency=8-30HZ')
       saveas(gca, 'RSA_permstat_tstat_HZ100_ngd_bp8_30.fig')
       saveas(gca, 'RSA_permstat_tstat_HZ100_ngd_bp8_30.jpg')
       
 
  %    imagesc(squeeze(p))   % p-map

%%

save(strcat(path_out,'_RSA_permtest_',num2str(win*1000),'ms_slide',num2str(slide*1000),'ms_time',num2str(t_start),'_',num2str(t_end),'wv_50HZ_permtest.mat'), 'neg_tsum', 'pos_tsum', 'h', 'p', 'c', 'stat', 'corr_all_allsubs', 'tval', 'tval_perm', 't1', 't2', 'p_pos', 'p_neg','mask_alpha', 'aaa', 'nrand', 'alpha', 'alpha_c')       