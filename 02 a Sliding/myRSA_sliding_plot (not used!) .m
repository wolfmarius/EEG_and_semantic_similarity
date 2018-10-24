% Only used for creating some plots to get an intuition of the data.



path_pre= 'F:\Masterarbeit_final\';

selvps = {'vp01'} %, 'vp02', 'vp03', 'vp04','vp05', 'vp06', 'vp07', 'vp08', 'vp09', 'vp10', 'vp11', 'vp12', 'vp13', 'vp14', 'vp15', 'vp16','vp17', 'vp18', 'vp19', 'vp20'};


for n=1:numel(selvps) 
   
   load(strcat(path_pre, selvps{n},'_slidRSA_steps200ms_slide20ms_time0_1HZ100.mat'))
        
   ind_on{n} =find(cat_sparse(:,1)== cat_sparse(:,2));
   ind_off{n}=find(cat_sparse(:,1)~= cat_sparse(:,2));
      
   auto_trial_all{n}  (:,:,:)= reshape([auto_trial{ind_on{n}}],numel(t1),numel(t1),[]);
   corr_all_allsubs{n}(:,:,:)= reshape([corr_all{:}],numel(t1),numel(t1),[]);  
   
end 

 t=(t1+t2)*0.5;  


%% Visualize

for n=1:numel(selvps)



tmp_ttest = corr_all_allsubs{n};
corr_all_allsubs_ttest =permute(tmp_ttest,[3,1,2]);

[h,p,c,stat]=ttest2(squeeze(corr_all_allsubs_ttest(ind_on{n},:,:)), squeeze(corr_all_allsubs_ttest(ind_off{n},:,:)) );

imagesc(t,t,squeeze(stat.tstat),[-3 3])
set(gca,'YDir','normal')
colorbar



end

%%



for n=1:numel(selvps)
    
figure('Name', strcat('figure_VP_nomeaning', num2str(n)),'NumberTitle','off');

subplot(2,2,1)

imagesc(t,t,squeeze(nanmean(corr_all_allsubs{n}(:,:,ind_on{n}),3)),[-0.05 0.05])
set(gca,'YDir','normal')
colorbar
title(strcat('ON','                              VP', num2str(n)))

subplot(2,2,2)
imagesc(t,t,squeeze(nanmean(corr_all_allsubs{n}(:,:,ind_on{n}),3)),[-0.05 0.05])
set(gca,'YDir','normal')
colorbar
title('OFF')

subplot(2,2,3)
diff_on_off=squeeze(nanmean(corr_all_allsubs{n}(:,:,ind_on{n}),3))- squeeze(nanmean(corr_all_allsubs{n}(:,:,ind_off{n}),3));
imagesc(t1,t1,diff_on_off,[-0.05 0.05])
set(gca,'YDir','normal')
colorbar
title('DIFF: ON-OFF')


tmp_ttest = corr_all_allsubs{n};
corr_all_allsubs_ttest =permute(tmp_ttest,[3,1,2]);

[h,p,c,stat]=ttest2(squeeze(corr_all_allsubs_ttest(ind_on{n},:,:)), squeeze(corr_all_allsubs_ttest(ind_off{n},:,:)) );
%ttest2, damit es independent ttets ist. Denn on und off sind independent

subplot(2,2,4)
imagesc(t,t,squeeze(stat.tstat),[-3 3])
set(gca,'YDir','normal')
colorbar
title('ttest')


%saveas(gca, strcat('figure_vp', num2str(n),'.jpg'))
%saveas(gca, strcat('figure_vp', num2str(n),'.fig'))

end

