function [corr_all, auto_trial]=vectorized_sliding_after_no_meaning(data,data_vec,c_ind,n_bins,cat_sparse, code_groups, code_blocks)     
                
             %select trials from categories
             c1=cat_sparse(c_ind,1);
             c2=cat_sparse(c_ind,2);          
                   
             trial_ind1=find(code_groups(:)==c1);
             trial_ind2=find(code_groups(:)==c2);                       
             
             trial_ind_tmp1 = [trial_ind1, code_blocks(trial_ind1)];
             trial_ind_tmp2 = [trial_ind2, code_blocks(trial_ind2)];             
  
                % vectorize correlations
                data_vec_tmp1=reshape(data_vec(trial_ind1,:,:),numel(trial_ind1)*n_bins,[]);
                data_vec_tmp2=reshape(data_vec(trial_ind2,:,:),numel(trial_ind2)*n_bins,[]);    
                            
               
                corr_tmp=corr(data_vec_tmp1',data_vec_tmp2','Type','Spearman');
                corr_tmp=0.5.*log((ones(size(corr_tmp))+corr_tmp)./(ones(size(corr_tmp))-corr_tmp));% fisher z-trans                             
                corr_tmp2=reshape(corr_tmp, numel(trial_ind1), n_bins, numel(trial_ind2), n_bins);
                                
              %  take care of same trial correlations
                        if c1==c2
                        on_tmp=eye(numel(trial_ind1)*n_bins);
                        on_tmp2=reshape(on_tmp,numel(trial_ind1), n_bins,numel(trial_ind1), n_bins);

                        same_tmp=zeros(size(on_tmp2));
                            for t_ind=1:numel(trial_ind1)
                                same_tmp(t_ind,:,t_ind,:)=1;
                            end
                        same_tmp2=same_tmp-on_tmp2;
                        rmsame=find(same_tmp);

                        takesame=find(ones(size(on_tmp2))-same_tmp2);

                        auto_trial_tmp=corr_tmp2;
                        auto_trial_tmp(takesame)=NaN;
                        auto_trial(:,:)=squeeze(nanmean(nanmean(auto_trial_tmp,1),3));

                        corr_tmp2(rmsame)=NaN;
                        else
                          auto_trial(:,:)=NaN;  
                            
                        end    
                        
               %  take care of same block correlations in different categories      
                     if c1~=c2   
                         for i=1:size(trial_ind1, 1)
                             for j=1:size(trial_ind2, 1)
                                 if trial_ind_tmp1(i,2) == trial_ind_tmp2(j,2)                                    
                                     %tmp1 = find(trial_ind_tmp1(:,1) == trial_ind_tmp1(i,1))                                     
                                     %tmp2 = find(trial_ind_tmp2(:,1) == trial_ind_tmp2(j,1))                             
                                     corr_tmp2(i,:,j,:)=NaN;
                                 end
                             end
                         end
                     end        
                      
                      corr_all(:,:)=squeeze(nanmean(nanmean(corr_tmp2,1),3));
     
end
  
         
      