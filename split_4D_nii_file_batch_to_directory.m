
filelist1=dir('**/*.nii');
list1={filelist1.('name')}';
pathlist={filelist1.('folder')}';
kmax=size(list1,1);



parfor k=1:kmax
    
strfile=list1{k,1};   
strfilepath=pathlist{k,1};   
totalfile=strcat(strfilepath,'\',strfile);
underlineLocations = strfind(strfile, '.');
% strfile2 = strfile(1:underlineLocations(4)-1);
% subject_name= strfile(1:underlineLocations(2)-1);
% condition_name=strfile(underlineLocations(3)+1:underlineLocations(4)-1);

strfile2 = strfile(1:underlineLocations(2)-1);
mkdir(strfile2)
output_directory_path=[strfilepath '\' strfile2]
spm_file_split(totalfile,output_directory_path );
end
