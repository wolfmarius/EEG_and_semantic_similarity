%% sliding RSA stats

addpath D:\Matlab_tools\fieldtrip-20170517
ft_defaults

%%


path_pre= 'C:\data\marius\02_RESULTS_sliding\';
path_out= 'C:\data\marius\03_RESULTS_permstat\';

selvps = {'vp01', 'vp02', 'vp03', 'vp04','vp05', 'vp06', 'vp07', 'vp08', 'vp09', 'vp10', 'vp11', 'vp12', 'vp13', 'vp14', 'vp15', 'vp16','vp17', 'vp18', 'vp19', 'vp20'};

nrand=1000;
alpha=0.05;
alpha_c=0.05;%p value for first level

%define matrix ind
% reorganize corr: off diagonal is symmetric, only use half of the matrix
for n=1:numel(selvps) 
    
 
 aaa= '_slidRSA_steps200ms_slide20ms_time0_1HZ100_hp30.mat'
 load(strcat(path_pre, selvps{n},aaa))
 
 
  clear all_cat
  [all_cat(:,1),all_cat(:,2)]=find(ones(n_groups));
        
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
       corr_all_allsubs{n}(:,:,i)=corr_tmp;     
              
   end
      
    n_group{n} = n_groups;
    
    ind_on{n}=find(all_cat(:,1)==all_cat(:,2));

    ind_off{n}=1:size(all_cat,1);
    ind_off{n}(ind_on{n})=[];
end



%% do randomisations
xx = 1;
for r=1:nrand   
     
    
 for n=1:numel(selvps)        
     
    ind_rand{n}= randperm(n_group{n}*n_group{n});  

    rand_on{n}  = ind_rand{n}(1:n_group{n});   
    rand_off{n} = ind_rand{n}(n_group{n}+1:end); 
    
    on_diag (n,:,:,:)= squeeze(nanmean(corr_all_allsubs{n}(:,:,rand_on{n}), 3));                 
    off_diag(n,:,:,:)= squeeze(nanmean(corr_all_allsubs{n}(:,:,rand_off{n}),3));
         
 end
 
   [h,p,c,stat]=ttest(on_diag, off_diag,'Alpha',alpha_c);  
 
     
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
        post=0  
    end
        
    pos_tsum{r}=sort(post,'descend');
    neg_tsum{r}=sort(negt,'ascend');
    
    %clear post negt      test mit und ohne
      
end
 

%%

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


% data
for n=1:numel(selvps)
    
    on_diag_notrand (n,:,:,:)= squeeze(nanmean(corr_all_allsubs{n}(:,:,ind_on{n}) ,3));                 
    off_diag_notrand(n,:,:,:)= squeeze(nanmean(corr_all_allsubs{n}(:,:,ind_off{n}),3));
          
end    

    [h,p,c,stat]=ttest(on_diag_notrand, off_diag_notrand ,'Alpha',alpha_c); % should I use ttest2 ??
    
 
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
     
%         while sig==1 i<min_pos
%        i=i+1;     
%        sig= data_post(i)>=pos_dist(i,round(alpha*nrand*(1/i)));
%        
%        if sig==1
%        sig_poscluster(i)=i;
%        mask_pos=mask_pos+(data_L_pos==ind_post(i));
%        end
% end
%         
        
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

% while sig==1 & i<min_neg
%        i=i+1;     
%        sig= data_negt(i)<=neg_dist(i,round(alpha*nrand*(1/i)));
%       if sig==1
%       sig_negcluster(i)=i;
%       mask_neg=mask_neg+(data_L_neg==ind_negt(i));
%       end
% end   



       %%
       % plot results
%         figure(10)
%         imagesc(mask_pos)
%         figure(11)
%         imagesc(mask_neg)
        
        mask_alpha=mask_pos+mask_neg;
        ind_z=find(mask_alpha==0);    
        mask_alpha(ind_z)=0.5;
       
     %%   
        figure()
   H= imagesc((t1+t2)*0.5,(t1+t2)*0.5,squeeze(stat.tstat),[-3 3]); 
       set(gca,'YDir','normal')
       colorbar
       set(H,'AlphaData',mask_alpha)
     %  title('meaningless,  sampling rate=100HZ,  frequency=1-7HZ')

%%
       
        figure()
   H= imagesc((t1+t2)*0.5,(t1+t2)*0.5, squeeze(nanmean(on_diag_notrand(:,:,:,:),1)) - squeeze(nanmean(off_diag_notrand(:,:,:,:),1)) ,[-0.05 0.05]);   
        set(gca,'YDir','normal')
        colorbar
        set(H,'AlphaData',mask_alpha)
        saveas(gca, 'permstat_on_minus_off_lp7.fig')
     
      
      
%%
      
save(strcat(path_out,'_permtest_',num2str(win*1000),'ms_slide',num2str(slide*1000),'ms_time',num2str(t_start),'_',num2str(t_end),'HZ100_hp30.mat'), 'neg_tsum', 'pos_tsum', 'h', 'p', 'c', 'stat', 'corr_all_allsubs', 'on_diag_notrand', 'off_diag_notrand', 'mask_alpha', 't1', 't2', 'aaa', 'nrand', 'alpha', 'alpha_c', 'p_neg', 'p_pos')       

      
      
      
      
      
      
      
      
      
      
      
      %% get peak of during stim with after stim correlation
%       
%      m1=nearest(t1,1);
%      m2=numel(t1);
%      
%      s1=nearest(t1,0);
%      s2=nearest(t1,0.5);
%      
%      figure(3)
%   imagesc(t1(m1:m2),t1(s1:s2), squeeze(stat.tstat(1,s1:s2,m1:m2)));
%      set(gca,'YDir','normal')
%     
%     figure(4)
%   plot(t1(s1:s2),  squeeze(nanmean(stat.tstat(1,s1:s2,m1:m2),3)));
%     
       %% get peak on diagonal (ti-ti), and get timecourse of tpeak with all others
%       
%       equal_t=eye(numel(t1));
%      equal_tind=find(equal_t);
%       tstat=squeeze(stat.tstat);
%       tstat_equal_t=tstat(equal_tind);
% 
%       figure(4)
%    plot(t1,tstat_equal_t);
%       
%       
%      ondiago=  on_diag_notrand;
%      ondiago_equal_t=ondiago(equal_tind);
%       
%      offdiago=  off_diag_notrand;
%      offdiago_equal_t=offdiago(equal_tind);
%     
%      figure(5)
%    plot(t1,ondiago_equal_t,t1,offdiago_equal_t)
%      
%      
%      [m,ind_max]=max(ondiago_equal_t-offdiago_equal_t);
%      
%      %ind_max=44     
%      for n=1:numel(selvps)
%          on_diag_tmp (n,:,:,:)=  squeeze(nanmean(corr_all_allsubs{n}(:, ind_max, ind_on{n}),3));    
%          off_diag_tmp (n,:,:,:)= squeeze(nanmean(corr_all_allsubs{n}(:, ind_max,ind_on{n}),3)); 
%      end
%      
%  %   max_on=squeeze(nanmean(nanmean(corr_all_allsubs(:,ind_max,:,ind_on),4),1));
%  %   max_off=squeeze(nanmean(nanmean(corr_all_allsubs(:,ind_max,:,ind_off),4),1));
%      
%      max_on = squeeze(nanmean(on_diag_tmp,1));
%      max_off= squeeze(nanmean(off_diag_tmp,1));
%      
%      figure(6)
%    plot(t1,max_on,t1, max_off)
%    
%    
%    if ind_max==1
%        ind_max=2;
%    end
%    
%      figure(7)
%    plot(t1,tstat(ind_max-1,:))
%    %plot(t1,tstat(43,:))