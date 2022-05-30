



function [meanIq,centrotype_vectors] =eCAPS_Tomas_kmeans_ICASSO(data, k_min, k_max, nbootstraps, data_ratio, kmeansreps)


EEGdata_one=data';  %% data = voxels x timepoints
M=nbootstraps; %number of boostraps
percent=data_ratio; %%proportion of resampled points e.g. 0.5 = 50%
reps=kmeansreps; %% number of kmeans repetitions





for k_number=k_min:k_max %% number of clusters

 



clear tic
clear resampled_data
clear k_model
clear concatenated_k_model_one
clear distance_matrix
clear initialProjection
clear projection_coordinates
clear dummy
clear linkage_clustering
clear sR_TOTAL
clear sR_model_one
clear sR_model_two
clear similarity_matrix_TOTAL
clear similarity_matrix_model_one 
clear distance_matrix_TOTAL
clear distance_matrix_model_two
clear distance_matrix_model_one
clear linkage_clustering_index_model_one 
clear linkage_clustering_index_TOTAL

    warning('off')

 %%rerun kmeans on randomly resampled data at 80%
 concatenated_k_model_one=[];
%  for i=1:M
    
for i=1:M
% parfor i=1:M %%use if you want parallel but large memory requirement
     tic
     
    warning('off')
    fprintf('%s\n',['Population 1 k=' num2str(k_number) ' : Resampling ' num2str(percent*100) '% of data Round ' num2str(i) '/' num2str(M) ]);  
 resampled_data_one=datasample(EEGdata_one,round(size(EEGdata_one,1)*percent),1);
%   resampled_data_one=EEGdata_one; %% no resampling

%%OPTION 1: kmeans clustering
    

%     disp(['kmeans on CPU for k=' num2str(k_number)]);
%       warning('off')
%  options = statset('UseParallel',0);
%    [~,k_model_one] = kmeans(resampled_data_one,k_number,'distance','correlation','replicates',reps,'empty','drop','MaxIter',1500, 'Options', options);
%     
%    
% %    %OR GPU based
  disp(['kmeans on GPU for k=' num2str(k_number)]);
  options = statset('UseParallel',0);
   [idx,k_model_one] = kmeans(gpuArray(resampled_data_one),k_number,'distance','correlation','replicates',reps,'empty','drop','maxiter',1500, 'Options', options);
    k_model_one=gather(k_model_one);
     idx=gather(idx);
     
%convert centroids to scale of data   
% for l=1:k_number    
% k_model_one(l,:) = mean(resampled_data_one(idx==l, :),1);
% end


%%%% OR preprocessing with PCA DIMENSIONALITY REDUCTION (much faster, but components do not seem that reliable)
%     warning('off')
%     disp(['PCA-kmeans on CPU for k=' num2str(k_number)]);
%     [pc,eigvec] = pca(resampled_data_one','NumComponents',k_number);
%     [idx,CP_reduced] = kmeans(pc,k_number,'distance','correlation','replicates',reps,'empty','drop','maxiter',1500);
% % %   convert centroids to scale of data  
% concatenated_CP=[];
% for l=1:k_number      
% CP = mean(resampled_data_one(idx==l, :));
% concatenated_CP=vertcat(concatenated_CP,CP);
% end
% k_model_one=concatenated_CP;
%     
    
    %% alternative PCA method but not as good in terms of cluster qulity: back projection (i.e. a different way of calculating the MS maps)
%     k_model_one = eigvec*CP_reduced; 
%      k_model_one=k_model_one';  



%%%OPTION 2: Non-negative matrix factorization (NNMF)
% disp('matlab NNMF on GPU');
% warning('off')
% 
% % % Matlab native NMF
% opt = statset('UseParallel',0, 'MaxIter',1000,'Display','off');
% [time_course, k_model_one]=nnmf(gpuArray(resampled_data_one),k_number, 'replicates',reps,'algorithm','als','options', opt);
%   k_model_one=gather(k_model_one);
  %%CPU only
%   [time_course, k_model_one]=nnmf(resampled_data_one,k_number, 'replicates',reps,'algorithm','als','options', opt);
  



% %% semi NMF version for negative (i.e. zscored) data (not great)
% [idx,Xout,Aout,Yout,numIter,tElapsed,finalResidual]=NMFCluster(resampled_data_one',k_number);
% k_model_one=Aout';


%%%OPTION 3: temporal ICA (tICA)
% disp('infomax temporal ICA');
% warning('off')
% [W,sphere] = runica(resampled_data_one', 'pca', k_number,  'extended', 0,'verbose' , 'off', 'lrate', 0.001);
% k_model_one = W * sphere;


% %%picard ICA instead? (faster)
%    [tmpdata,eigvec] = pca(resampled_data_one','NumComponents',k_number);
%      [~, W_] = picard(tmpdata', 'verbose', false, 'mode', 'standard', 'tol', 1e-5);
%       k_model_one = W_*pinv(eigvec);

          
% %%%OPTION 4: spatial ICA (sICA)
% disp('infomax spatial ICA');
% warning('off')
% [W,sphere] = runica(resampled_data_one, 'pca', k_number,  'extended', 0,'verbose' , 'off', 'lrate', 0.001);
% k_model_one = (resampled_data_one'*pinv(W * sphere))';  %% data*inv(activations)

% % % %%picard ICA instead? (faster)
%    [~,eigvec] = pca(resampled_data_one','NumComponents',k_number);
%     [~, W_] = picard(eigvec', 'verbose', false, 'mode', 'standard', 'tol', 1e-5);
%     k_model_one=W_*eigvec';
% 
% 



 toc 
concatenated_k_model_one=vertcat(concatenated_k_model_one,k_model_one); 
end

 


   
  concatenated_k_model_one_double=concatenated_k_model_one;


 
 % % MODEL_ONE calculate centrotypes
distance_matrix_model_one=squareform(pdist(concatenated_k_model_one_double,@Tomas_microstatetopo_sqrtdistance));
distance_matrix_model_one_vectorized=squareform(distance_matrix_model_one,'tovector');
similarity_matrix_model_one=squareform(pdist(concatenated_k_model_one_double,@Tomas_microstatetopo_adjacency));
[sR_model_one.cluster.partition,sR_model_one.cluster.dendrogram.Z]=hcluster_Tomas(distance_matrix_model_one,'ward');
sR_model_one.cluster.similarity=similarity_matrix_model_one;
centrotype_index_model_one=icassoIdx2Centrotype(sR_model_one,'partition',sR_model_one.cluster.partition(k_number,:)); 
centrotype_vectors=concatenated_k_model_one_double(centrotype_index_model_one, :);

for i=1:30
 linkage_clustering_index_model_one(:,i)= cluster(sR_model_one.cluster.dendrogram.Z,'maxclust',i);
end

%%ICASSO


[Iq2,in_avg,ext_avg,in_min,ext_max] =icassoStability(sR_model_one,k_number,'none');
% Iq2
meanIq(:,k_number)=mean(Iq2);


%%Silhouette score


% linkage_eva_one = evalclusters(concatenated_k_model_one,linkage_clustering_index_model_one,'silhouette','ClusterPriors', 'empirical','Distance', distance_matrix_model_one_vectorized);
% % figure;plot(linkage_eva_one.InspectedK,linkage_eva_one.CriterionValues);title('linkage');
% fprintf ('\n');
% disp( ['k=' num2str(k_number) '  ' title1 ', Silhouette Value = '   num2str(linkage_eva_one.CriterionValues(k_number))]);
% silhouette_value_model_one(:,k_number)=linkage_eva_one.CriterionValues(k_number);
% 
% linkage_eva_one.CriterionValues
% size(linkage_eva_one.CriterionValues)



% [ SOptimal ] = multislice_pair_labeling(horzcat(linkage_clustering_index_model_one(:,k_number),linkage_clustering_index_model_two(:,k_number)));
% SOptimal_model_one=SOptimal(:,1); SOptimal_model_two=SOptimal(:,2);

%%OR spectral clustering?

title1=['eCAPS cluster quality Iq= ' num2str(meanIq(k_number)) ' :  k=' num2str(k_number) ' data ratio=' num2str(data_ratio) ' resamples=' num2str(nbootstraps)];


%%CCA projection TOTAL
alpha=0.7;
CCAradius=3*max(std(distance_matrix_model_one)); %default
epochs=10;
outputDimension=2;
initialProjection_TOTAL=mmds(distance_matrix_model_one); initialProjection_TOTAL=initialProjection_TOTAL(:,1:2);  %MMDS initiliazation
dummy=rand(size(distance_matrix_model_one,1),outputDimension);
projection_coordinates_TOTAL=cca(dummy,initialProjection_TOTAL,epochs,distance_matrix_model_one,alpha,CCAradius);

% %CCA projection MODEL_ONE
figure; 
gscatter(projection_coordinates_TOTAL(1:(end),1),projection_coordinates_TOTAL(1:(end),2),linkage_clustering_index_model_one(:,k_number),[],'.' );title(title1);
hold on; gscatter(projection_coordinates_TOTAL(centrotype_index_model_one,1),projection_coordinates_TOTAL(centrotype_index_model_one,2),[],'k','.',25,'MarkerFaceColor','k');
lgd=legend(); lgd.String(k_number+1)={'centroids'};
xlim([min(projection_coordinates_TOTAL(:,1)) max(projection_coordinates_TOTAL(:,1))]); ylim([min(projection_coordinates_TOTAL(:,2)) max(projection_coordinates_TOTAL(:,2))]);

% %%CCA projection MODEL TOTAL
% group=ones(size(projection_coordinates_TOTAL,1),1);
% group(size(projection_coordinates_TOTAL,1)/2:end)=2;
% % 
% figure; gscatter(projection_coordinates_TOTAL(:,1),projection_coordinates_TOTAL(:,2),group,'k', 'op'); %%%legend_labels=[{'ec'}, {'eo'}];legend(legend_labels);
% hold on; gscatter(projection_coordinates_TOTAL(:,1),projection_coordinates_TOTAL(:,2),linkage_clustering_index_TOTAL(:,k_number));
% hold on; gscatter(projection_coordinates_TOTAL(centrotype_index_TOTAL,1),projection_coordinates_TOTAL(centrotype_index_TOTAL,2),[],'k','.',25,'MarkerFaceColor','k');
% lgd=legend(); lgd.String(k_number+1)={'centroids'};

% hold on; gscatter(projection_coordinates_TOTAL(centrotype_index_model_one,1),projection_coordinates_TOTAL(centrotype_index_model_one,2),[],'k','x',15,'MarkerFaceColor','k');
% hold on; gscatter(projection_coordinates_TOTAL(centrotype_index_model_two_in_TOTAL,1),projection_coordinates_TOTAL(centrotype_index_model_two_in_TOTAL,2),[],'k','x',15,'MarkerFaceColor','k');
% 
% 


%%plot microstate topographies TOTAL
% load('nbt_coolWarm');
% colormap(coolWarm)
%  
% figure;
% for i=1:k_number
%     
%     hold(subplot(2,k_number,i),'on'); title(i)
%     topoplot(centrotype_vectors_model_one(i,:),ALLEEG(1).chanlocs,'plotrad',0.65,'headrad','rim','electrodes','off', 'colormap', coolWarm,'numcontour',0);
%      
% 
%     
% 
% end
% 

end


% figure; hb = bar(1:1:length(meanIq), meanIq(:), 'grouped');
% title(title1);

