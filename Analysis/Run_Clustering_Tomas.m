function [CP2,Disp,Std_Clusters,idx,d,sfrac] = Run_Clustering_Tomas(XONn,n_clusters,mask,brain_info,maskP,maskN,n_rep,idx_sep_seeds,SeedType)

    fprintf('\n')
    disp('Clustering data...please wait')
    % Number of seeds
    n_seeds = size(idx_sep_seeds,3);
    
    % Number of subjects
    n_subjects = size(idx_sep_seeds,2);

    % Number of possible seed combinations
    switch n_seeds
        case 1
            n_combos = 1;
            combvec = [1];
        case 2
            n_combos = 3;
            combvec = [1 0; 0 1; 1 1];
        case 3
            n_combos = 7;
            combvec = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1];
        otherwise
            errordlg('PROBLEM AT SEED FRACTIONS');
    end

    % Will contain the fraction of frames linked to a given seed (if using
    % the 'Intersection' method)
    sfrac = zeros(n_subjects,n_clusters,n_combos);
    
    
  
    
    % 'Filtering so that we only use the largest activation and deactivation
    % spots for clustering
%     XONn_filtered = CAP_mask4kmeans(XONn,maskP,maskN,6,mask,brain_info);

%     XONn_filtered=XONn; %%Tomas skip masking
%     max(XONn_filtered(:))
%     min(XONn_filtered(:))
    
%     size(XONn_filtered)  %% XONn_filtered == voxels x timepoints
    
% %%%low-pass filtering (optional)
% fpass=1 ; %% upper cutoff in Hz
% fs=40; %%sampling rate
% XONn_filtered=(lowpass(XONn_filtered',fpass,fs))';
    
% 
%      XONn_filtered=zscore(XONn, 0, 2); %%Tomas zscore ?
%      XONn=zscore(XONn, 0, 2); %%Tomas zscore

     
% size(XONn_filtered')
% max((isnan(XONn_filtered(:)')))
    
    % Rows datapoints, columns variable (so here every 70700 activation is
    % a datapoint)
    % idx will contain 1462 elements (the index of the cluster to which the
    % considered datapoint belongs
    
    
    

    
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  kmeans on CPU  %%%%%%%%%%%%%%%%
% disp('kmeans on CPU');  %%N.B. kmeans provides negative CAPs maps
% tic %%Tomas edit
%     options = statset('UseParallel',0);
%     [idx,CP] = kmeans(XONn',n_clusters,'distance','correlation','replicates',n_rep,'empty','drop','maxiter',1000, 'Options',options);
% 
% %convert centroids to scale of data   
% for l=1:max(idx)       
% CP(l,:) = mean(XONn(:,idx==l),2);
% end
%     
% for time=1:length(idx)
% % time_course(:,time)=nnls(CP', XONn(:,time));
% time_course(:,time)=tntnn(CP', XONn(:,time));
% end
% toc
% 
% max(CP')
% min(CP')
% 
% max(time_course(:))
% min(time_course(:))   


    %%GPU version
    disp('kmeans on GPU');
    tic %%Tomas edit
    options = statset('UseParallel',0);  %%out of memory if run in parallel mode
    XONn=gpuArray(XONn);
% %      class(XONn)
    [idx,CP] = kmeans(XONn',n_clusters,'distance','correlation','replicates',n_rep,'empty','drop','maxiter',500, 'Options',options);
    idx=gather(idx);
    CP=gather(CP);
    XONn=gather(XONn);
    toc

    
    
%convert centroids to scale of data   
for l=1:max(idx)       
CP(l,:) = mean(XONn(:,idx==l),2);
end
%     
% for time=1:length(idx)
% % time_course(:,time)=nnls(CP', XONn(:,time));
% time_course(:,time)=tntnn(CP', XONn(:,time));
% end

%     

%%%%%%%%%%%%%%%%% Matlab's non-negative matrix factorisation  %%%%%%%%%%%%%
% disp('matlab NNMF on CPU');
% tic
% warning('off')
% opt = statset('UseParallel',1, 'MaxIter',1000,'Display','final');
% [time_course, CP]=nnmf(XONn',n_clusters, 'replicates',6,'algorithm','als','options', opt);
% % [CP,time_course]=CauchyNMF(XONn,n_clusters,'SCALE',-2,'VERBOSE',1);time_course=time_course';CP=CP';
% % max(time_course)
% % min(time_course)
% time_course=time_course';
% % size(time_course)
% % size(CP)
% [~,idx]=max(time_course,[],1);
% idx=idx';
% toc


% %%NNMF CLUSTERING
% disp('NMF clustering on CPU');
% warning('off')
% tic
% % option.algorithm='vsmf';
% option.algorithm='nmfrule';
% % option.algorithm='nmfnnls';
% % option.algorithm='sparsenmfnnls';
% % option.algorithm='convexnmfrule';
% % option.algorithm='orthnmfrule';
% % option.algorithm='sparsenmf2rule';
% % option.algorithm='sparsenmf2nnqp';
% 
% % option.t1=1;  %% A is non-negative i.e. CAPs maps
% % option.t2=1;  %% Y is non-negative i.e. time-series
% % % option.reorder=0;
% 
% % max(XONn(:))
% % min(XONn(:))
% [idx,Xout,Aout,Yout,numIter,tElapsed,finalResidual]=NMFCluster(XONn,n_clusters, option);
% % [idx,Xout,Aout,Yout,numIter,tElapsed,finalResidual]=NMFCluster(XONn,n_clusters);
% 
% % size(Xout)  %% reconstructed XONn
% % size(Aout)  %% CAPs maps == voxels x components
% % size(Yout)  %% time-series == components x timepoints
% CP=Aout';
% time_course2=Yout;
% toc

% sampling_rate=40; %%in Hz
% time_in_seconds=(1:length(idx))./sampling_rate; 
% figure;plot(time_in_seconds(1:1600)',time_course2(1:2, 1:1600)'); %%1600 x @40 Hz = 40 seconds
% component_correlations2=corr(time_course2');
% figure;heatmap(component_correlations2, 'Colormap',flipud(cbrewer('div','RdBu',10)), 'ColorLimits',[-1 1]);

% for time=1:length(idx)
% % time_course(:,time)=nnls(CP', XONn(:,time));
% time_course(:,time)=tntnn(CP', XONn(:,time));
% end



%% TEMPORAL SMOOTHING (to be consistent with EEG microstate analysis)
%% option1: median filtering
%   idx=medfilt1( idx, 3 ); %%Matlab native median filter
%   idx=medfilt1m( idx, 1 ); %%adaptive median filter
%   idx=idx';

% % option2: smoothing via Koenig microstate toolbox
% opts.b=2; %%smoothing window in samples (@40Hz 1 sample=25ms)
% opts.lambda=5; %smoothing weight
% idx = smoothing(XONn,CP',n_clusters,opts);
% idx=idx';


%% option3: reject small segments via Poulsen microstate toolbox
% opts.minTime = 1; %%minimum segment size in samples
% idx = MicroSmooth(XONn, CP', 'reject segments' ,opts); %% 'reject segments' (default) or 'windowed'.
% idx=idx';

% % option4: windowed smoothing via Poulsen microstate toolbox
% opts.b=2; %%smoothing window in samples (@40Hz 1 sample=25ms)
% opts.lambda=5; %smoothing weight
% idx = MicroSmooth(XONn, CP', 'windowed' ,opts); %% 'reject segments' (default) or 'windowed'.
% idx=idx';



    % idx2counts is of size K (number of clusters) and has the number of
    % datapoints classified within a given cluster)

    % disp('idx2counts:');
%     idx2counts = histc(idx, 1:max(idx));
%     
%     % Output = Input(IX)
%     [~,IX] = sort(idx2counts,'descend');
%     
%     % Size Kx70700 (location of each cluster); clusters are put with 'the
%     % most prominent one first'
%     CP = CP(IX,:); idx2 = idx; % order by occurrence
% 
%     % Changes the datapoint indexes so that they fit the new clusters order
%     for l=1:max(idx), idx2(idx==IX(l))=l; end
%     idx=idx2;

%     CP2 = zeros(n_clusters,size(CP,2));
    CP2=CP; %%Tomas no cluster sorting
    
    Disp = zeros(1,n_clusters);
    Std_Clusters = zeros(size(CP,2),n_clusters);

    % For each cluster index
    for l=1:max(idx)       
        % Averages all data points belonging to one specific cluster, and
        % stores the obtained pattern as a cell in CP2
%         CP2(l,:) = mean(XONn(:,idx==l),2); %./ ( std(XON(:,idx==l),[],2) ) * sqrt(length(idx==l)); % Liu&Duyn

        % Measure of dispersion within the cluster considered
        Disp(l) = mean(corr(CP2(l,:)',XONn(:,idx==l)));

        Std_Clusters(:,l) = std(XONn(:,idx==l),[],2);
    end
    
    % d contains the correlation values of all frames to the CAPs
    r = corr(CP2',XONn);

    d = zeros(1,length(idx));
    for k=1:max(idx)
        d(idx==k) = r(k,idx==k);
    end
    
    % Added part to compute the fraction of frames assigned to a given seed
    % if using the Union option (in which a data point is retained
    % as long as at least one seed region becomes significantly (de)active)
    if strcmp(SeedType,'Union')
        
        % idx_all will contain the clustering indices as put on a whole
        % temporal scale (time x subjects)
        idx_all = zeros(size(idx_sep_seeds,1),n_subjects);
        
        % Index to properly fill in the matrix by adding up the number of
        % frames per subject every time
        tmp_loc = 1;
        

        for s = 1:n_subjects
            tmp = sum(squeeze(idx_sep_seeds(:,s,:)),2);
            tmp(tmp >= 1) = 1;
            tmp = logical(tmp);
            idx_all(tmp,s) = idx(tmp_loc:(tmp_loc+sum(tmp)-1));
            tmp_loc = tmp_loc + sum(tmp);
        end
            
        
        
        % I will compute my fractions for each seed combination of
        % interest; for example, 3 seeds would yield 7 possible
        % combinations
        for s = 1:n_subjects
            for t = 1:size(idx_sep_seeds,1)
                
                tmp = squeeze(idx_sep_seeds(t,s,:));   
                
                % If there is at least one seed active or deactive at
                % the point of interest, we update the sfrac count at the
                % appropriate CAP
                if sum(tmp) > 0
                    sfrac(s,idx_all(t,s),find(ismember(combvec,tmp','rows'))) = ...
                        sfrac(s,idx_all(t,s),find(ismember(combvec,tmp','rows'))) + 1;   
                end
            end
        end 
    end 
    
%%Tomas compute temporal correlation between components
% time_course = time_course(IX,:); %% sort time course above

% sampling_rate=40; %%in Hz
% time_in_seconds=(1:length(idx))./sampling_rate; 
% figure;plot(time_in_seconds(1:1600)',time_course(1:2, 1:1600)'); %%1600 x @40 Hz = 40 seconds
% component_correlations=corr(time_course');
% figure;heatmap(component_correlations, 'Colormap',flipud(cbrewer('div','RdBu',10)), 'ColorLimits',[-1 1]);

%%%TOMAS EDIT %% variance explained by each component
% for comp=1:n_clusters
%     data_variance = sum(var(XONn));
%     projected_variance(comp) = sum(var(CP2(comp,:)' * time_course(comp, :)));
%     variance_explained (comp)= 100* projected_variance(comp)/data_variance;
% end  
% variance_explained=double(variance_explained)
% global_variance_explained=sum(variance_explained)

end