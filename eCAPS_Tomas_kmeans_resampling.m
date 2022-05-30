%%% CCA directly on source EEG data (e.g. voxels x timepoints)
%loasd data
% tic




function silhouette_value_model_one=eCAPS_Tomas_kmeans_resampling(data, k_min, k_max, nbootstraps, data_ratio, kmeansreps)


EEGdata_one=data';  %% data = voxels x timepoints
M=nbootstraps; %number of boostraps
percent=data_ratio; %%proportion of resampled points e.g. 0.5 = 50%
reps=kmeansreps; %% number of kmeans repetitions


title1=['eCAPS kmeans reps=' num2str(reps)];

for k_number=k_min:k_max %% number of clusters

  
 



clear tic
clear resampled_data
clear k_model
clear concatenated_k_model_one
clear concatenated_k_model_two
clear distance_matrix
clear initialProjection
clear projection_coordinates
clear dummy
clear linkage_clustering
clear sR_TOTAL
clear sR_model_one
clear sR_model_two
clear similarity_matrix_TOTAL
clear similarity_matrix_model_two
clear similarity_matrix_model_one 
clear distance_matrix_TOTAL
clear distance_matrix_model_two
clear distance_matrix_model_one
clear linkage_clustering_index_model_one 
clear linkage_clustering_index_model_two 
clear linkage_clustering_index_TOTAL

 %%rerun kmeans on randomly resampled data at 80%
 concatenated_k_model_one=[];
 for i=1:M
    fprintf('%s\n',['Group 1: Randomization using kmeans Round ' num2str(i) '/' num2str(M) ]);  
 resampled_data_one=datasample(EEGdata_one,round(size(EEGdata_one,1)*percent),1);
%   resampled_data_one=EEGdata_one; %% no resampling

%%kmeans clustering
    disp(['kmeans on GPU for k=' num2str(k_number)]);
    tic %%Tomas edit
    [~,k_model_one] = kmeans(gpuArray(resampled_data_one),k_number,'distance','correlation','replicates',reps,'empty','drop','maxiter',500);
%       [~,k_model_one] = kmeans(resampled_data_one),k_number,'distance','correlation','replicates',reps,'empty','drop','maxiter',500);
%     size(k_model_one)
    k_model_one=gather(k_model_one);
    toc
    
    
  
 concatenated_k_model_one=vertcat(concatenated_k_model_one,k_model_one);

 end
  
 


   
  concatenated_k_model_one_double=concatenated_k_model_one;


 
 % % MODEL_ONE calculate centrotypes
distance_matrix_model_one=squareform(pdist(concatenated_k_model_one_double,@Tomas_ecaps_sqrtdistance));
distance_matrix_model_one_vectorized=squareform(distance_matrix_model_one,'tovector');
similarity_matrix_model_one=squareform(pdist(concatenated_k_model_one_double,@Tomas_ecaps_adjacency));
[sR_model_one.cluster.partition,sR_model_one.cluster.dendrogram.Z]=hcluster_Tomas(distance_matrix_model_one,'ward');
sR_model_one.cluster.similarity=similarity_matrix_model_one;
centrotype_index_model_one=icassoIdx2Centrotype(sR_model_one,'partition',sR_model_one.cluster.partition(k_number,:)); 
centrotype_vectors_model_one=concatenated_k_model_one_double(centrotype_index_model_one, :);
for i=1:12
 linkage_clustering_index_model_one(:,i)= cluster(sR_model_one.cluster.dendrogram.Z,'maxclust',i);
end


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


figure; hb = bar(1:1:length(silhouette_value_model_one), silhouette_value_model_one(:), 'grouped');
legend(title1);

