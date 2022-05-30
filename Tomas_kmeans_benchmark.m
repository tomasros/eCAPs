%%tomas clustering benchmark

k=5;
replicates=5;
testMatrix=randn(10000);
tic
fprintf('\n')
disp('serial CPU')
testcluster=kmeans(testMatrix,k, 'replicates',replicates);
toc

tic
fprintf('\n')
disp('parallel CPU')
options = statset('UseParallel',1);
testcluster2=kmeans(testMatrix,k,'replicates',replicates, 'Options',options);
toc


tic
fprintf('\n')
disp('GPU')
testcluster3=kmeans(gpuArray(testMatrix),k,'replicates',replicates);
testcluster3=gather(testcluster3);
toc

tic
fprintf('\n')
disp('parallel CPU + GPU')
options = statset('UseParallel',1);
testcluster4=kmeans(gpuArray(testMatrix),k,'replicates',replicates, 'Options',options);
toc