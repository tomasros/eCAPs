%generate_masks_from_atlas.m
%This script generate mask files from any atlases you prefer.
%K. Nemoto 25 April 2015


xA=spm_atlas('load');
S=spm_atlas('select',xA);

for i = 1:size(S,2)
    fname=strcat(S{i},'.nii');
    VM=spm_atlas('mask',xA,S{i});
    VM.fname=fname;
    spm_write_vol(VM,spm_read_vols(VM));
end
