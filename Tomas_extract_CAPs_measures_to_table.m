%% extract MS measures from CAPs toolbox output

%% first load CAPS output *.mat file into workspace
%% info is in Outputs variable
disp('Saving eCAPs data....')
clear eCAPs
clear BIOMARKER_table
clear CLASSIFICATION_matrix

project_name = OverallInfo.ProjectTitle;
directory=OverallInfo.SaveDir;


k_number=Parameters.KMeansClustering.NumberClusters; %%total cluster number
group_number=length(OverallInfo.NumberSubjects);



index=0;
for g=2:group_number
    group_subject_number=OverallInfo.NumberSubjects{1,g};
    
    for s=1:group_subject_number
        index=index+1;
        entry(index,1)=[g];
        entry(index,2)=[s];
        
      
    eCAPs(index).Group=entry(index,1);    
    eCAPs(index).Subjects=entry(index,2);
    
    
% eCAPs(index).SubjectID=OverallInfo.SubjectID{1,g}{s,1:k_number};

eCAPs(index).Time_coverage=Outputs.Metrics.Occurrences{1, g}.frac.state(s,1:k_number);
eCAPs(index).Duration=Outputs.Metrics.AverageExpressionDuration{1, g}(s, 3:2+k_number); 
eCAPs(index).Occurrence=Outputs.Metrics.NumberEntries{1, g}(s, 3:2+k_number);

    end
end
   
 BIOMARKER_table =struct2table(eCAPs);

 
toolbox_file_location=which('CAP_TB.m');
toolbox_path=toolbox_file_location(1:end-8);
cd(toolbox_path) %%change directory
table_name = OverallInfo.ProjectTitle;
s = strcat(toolbox_path, 'SavedData\', table_name, '.csv');
 writetable(BIOMARKER_table,s); 

if group_number>1 && group_subject_number>1
    
    
    savetext = fopen([OverallInfo.SaveDir '\'  project_name '_effect_size.txt'],'w');
    
    
    
%%plot Time Coverage results
for K=1:k_number
% figure(100+K)
[p,tbl,anovastats]  = anova1(BIOMARKER_table.Time_coverage(:,K),BIOMARKER_table.Group, 'off');
% [c,~,~,gnames] = multcompare(anovastats);
[th,tp,tci,tstats] = ttest2(BIOMARKER_table.Time_coverage(BIOMARKER_table.Group==2,K),BIOMARKER_table.Time_coverage(BIOMARKER_table.Group==3,K));
cohens_d_effect_size=(anovastats.means(2)-anovastats.means(1))/tstats.sd;

% title(['Time Coverage for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant'))  ]);
disp(['Time Coverage for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant'))  ]);
fprintf(savetext, ['Time Coverage for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant')) '\n' ]);
end


% % plot Duration results
for K=1:k_number
% figure(200+K)
[p,tbl,anovastats]  = anova1(BIOMARKER_table.Duration(:,K),BIOMARKER_table.Group, 'off');
% [c,~,~,gnames] = multcompare(anovastats);
[th,tp,tci,tstats] = ttest2(BIOMARKER_table.Duration(BIOMARKER_table.Group==2,K),BIOMARKER_table.Duration(BIOMARKER_table.Group==3,K));
cohens_d_effect_size=(anovastats.means(2)-anovastats.means(1))/tstats.sd;

% title(['Duration for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant'))  ]);
disp(['Duration for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant'))  ]);
fprintf(savetext, ['Duration for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant')) '\n' ]);
end

% plot Occurrence results
for K=1:k_number
% figure(300+K)
[p,tbl,anovastats]  = anova1(BIOMARKER_table.Occurrence(:,K),BIOMARKER_table.Group, 'off');
% [c,~,~,gnames] = multcompare(anovastats);
[th,tp,tci,tstats] = ttest2(BIOMARKER_table.Occurrence(BIOMARKER_table.Group==2,K),BIOMARKER_table.Occurrence(BIOMARKER_table.Group==3,K));
cohens_d_effect_size=(anovastats.means(2)-anovastats.means(1))/tstats.sd;

% title(['Occurrence for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant'))  ]);
disp(['Occurrence for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant'))  ]);
fprintf(savetext, ['Occurrence for eCAP ' num2str(K) '  :  Effect size= ' num2str(round(cohens_d_effect_size,1, 'significant')) char(9) '   p= ' num2str(round(tp,1, 'significant')) '\n' ]); 
end

else    
    
end

    
    %%%%%%%%%%%%%%%%%%SAVE CAPs maps to images

figure
parfor i =1:k_number
map_filename=[OverallInfo.SaveDir '\' project_name '_eCAP' num2str(i) '.nii'];
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',map_filename,'options4Brainnetviewer.mat');
title(['eCAP ' num2str(i)]);
frame=getframe(gcf);
 im{i} = frame2im(frame);
 new=frame.cdata;
imwrite(new,[ map_filename '.tif'])
end

%%save as gif
filename = [pwd '\SavedData\' project_name '_animated_eCAPs.gif']; % Specify the output file name
for i = 1:k_number
    [A,map] = rgb2ind(im{i},256, 'nodither'); 
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.8);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.8);
    end
   
end

 BIOMARKER_table=splitvars(BIOMARKER_table);
 assignin('base','eCAPs', eCAPs);
 assignin('base','BIOMARKER_table', BIOMARKER_table);
 
 predictorNames=BIOMARKER_table.Properties.VariableNames(3:end);
 CLASSIFICATION_cell = table2cell(BIOMARKER_table(:, predictorNames));
 CLASSIFICATION_matrix=cell2mat(CLASSIFICATION_cell);
 assignin('base','CLASSIFICATION_matrix', CLASSIFICATION_matrix);
 
group1_name=2;
assignin('base','group1_name', group1_name);
 
group2_name=3;
assignin('base','group2_name', group2_name);
 
 

% fclose(savetext)
close   

