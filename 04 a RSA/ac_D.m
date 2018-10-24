function [ac_D, cs_D_names] = ac_D(uniq_words, D)

%used in script RSA_preprocess. It loads the RDMs and brings it into
%all-cat form (as opposed to cat-sparse form)



if exist('D_wv')
    D = D_wv;
end

%% Bring unique words and D in cat-sparse form

   n_words = numel(uniq_words);
   corr_mat_cat=zeros(n_words);
   
for ind=1:n_words       %this loop creates a square matrix with 1 on the diagonal and in the upper triangle
           corr_mat_cat(:,ind)=[ones(ind,1);zeros(n_words-ind,1)];
end

         [cat1,cat2]=find(corr_mat_cat); %returns row and column index of non-zero elements
          cat_sparse_D=[cat1,cat2];

          
% unique words
for c_ind=1:size(cat_sparse_D,1)
     c1=cat_sparse_D(c_ind,1);
     c2=cat_sparse_D(c_ind,2); 
     
     cs_D_names{c_ind} = [uniq_words{c1}, uniq_words{c2}];      
          
end           
          
    
%% strecht out D, ignoring NaN   

for i=1:size(D,1)     %if the input D is NaN on the Diagonal. Otherwise you can take it out
    D(i,i)= 0;
end

%%
idx = tril(true(size(D)),-1);  % take the -1 out if you want as well the diagonal to be NaN
D(idx) = NaN;

D_stretch = D(:);       
D_stretch(~any(~isnan(D_stretch), 2),:)=[];
ac_D = D_stretch';

end

