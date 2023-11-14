function ME_spm8_skullstrip(infile,outroot,thresh,ndilate_erodes)

if (nargin < 1), infile  = ''; end
if (nargin < 2), outroot = ''; end
if (nargin < 3), thresh  = 50; end
if (nargin < 4), ndilate_erodes  = 2; end

if (isempty(infile))
    if (ispref(mfilename(),'niftipath')), niftipath = getpref(mfilename(),'niftipath');
    else, niftipath = [pwd() filesep()]; end
    [infile,niftipath] = uigetfile('*.nii','Select the NIFTI to skull strip',niftipath);
    if (~ischar(infile)), return; end                   % user hit cancel
    setpref(mfilename(),'niftipath',niftipath);         % remember for next time
    infile = [niftipath infile];
end

if (isempty(outroot))
    [~, outroot] = fileparts(infile);
end

% --- make input filename have a fullpath ---
[inpath,name,ext] = fileparts(infile);
if (isempty(inpath)), inpath = pwd(); end
opwd = pwd();
cd(inpath);     
inpath = pwd();
infile = [inpath filesep() name ext]; % this turns relative path into absolute path
cd(opwd);

% --- find the path to the tissue masks ---
spmfile = which('spm');
spmpath = fileparts(spmfile);
tissuepath = [spmpath filesep() 'tpm' filesep()];

% --- set up SPM batch params ---
matlabbatch{1}.spm.spatial.preproc.data = {[infile ',1']};
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.opts.tpm = {[tissuepath 'grey.nii'] [tissuepath 'white.nii'] [tissuepath 'csf.nii']};
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2 2 2 4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};

% --- execute segmentation of infile ---
inputs = cell(0, 1);
spm('defaults', 'FMRI');    % not sure if this is really needed
spm_jobman('initcfg');    % This will get auto-run by SPM, just once
res = spm_jobman('serial', matlabbatch, '', inputs{:});

% --- build skullstripped (and bias corrected) from segmented brain tissues ---
v0 = spm_vol(infile);               % input file
v1 = spm_vol(res{1}.c1{1});         % tissue classes
v2 = spm_vol(res{1}.c2{1});
v3 = spm_vol(res{1}.c3{1});
v4 = spm_vol(res{1}.biascorr{1});   % bias corrected

I0 = spm_read_vols(v0);
I1 = spm_read_vols(v1);
I2 = spm_read_vols(v2);
I3 = spm_read_vols(v3);
I4 = spm_read_vols(v4);

% --- make skull stripped mask ---
Imask  = (I1 + I2 +I3)>(thresh/100); % mask
Imask  = imfill(Imask,26,'holes');   % fill holes
% --- Further fill holes w/ dilates followed by erodes
for i=1:ndilate_erodes, Imask  = imdilate(Imask,ones(3,3,3)); end
for i=1:ndilate_erodes, Imask  = imerode(Imask,ones(3,3,3)); end

% --- make masked raw image and masked bias corrected image ---
Istrip = I0 .* Imask;                % skull stripped
Ibias  = I4 .* Imask;                % skull stripped & bias corrected
vx     = v0;                         % copy volume info (will overwrite image pixel data)

% --- Write NIFTIs ---
% filename1   = [inpath filesep() outroot '_SPM_mask.nii'];
% vx.fname   = filename1;
% vx.descrip = 'SPM Skull stripped mask';
% fprintf(1,'Writing NIFTI: %s\n',filename1);
% spm_write_vol(vx,Imask);

% filename2  = [inpath filesep() outroot '_SPM_masked.nii'];
% vx.fname   = filename2;
% vx.descrip = 'SPM Skull stripped';
% fprintf(1,'Writing NIFTI: %s\n',filename2);
% spm_write_vol(vx,Istrip);

filename3   = [inpath filesep() outroot '_SPM_biascor.nii'];
vx.fname   = filename3;
vx.descrip = 'SPM bias corrected';
fprintf(1,'Writing NIFTI: %s\n',filename3);
spm_write_vol(vx,I4);

% filename4 = [inpath filesep() outroot '_SPM_biascor_masked.nii'];
% vx.fname   = filename4;
% vx.descrip = 'SPM Skull stripped and bias corrected';
% fprintf(1,'Writing NIFTI: %s\n',filename4);
% spm_write_vol(vx,Ibias);

% --- Delete temp files ---
delete(res{1}.c1{1});
delete(res{1}.c2{1});
delete(res{1}.c3{1});
delete(res{1}.biascorr{1});
delete(res{1}.snfile{1});
delete(res{1}.isnfile{1});

end
