
function [vol,hdr,filename4D] = dicom_spm2nifti(infiles,outfile,make4D)
% Use SPM to convert Dicoms to Nifti volume
%
% "infiles" can be cell array of filenames, OR cell array of dicominfo() structs
%
% M.Elliott 9/13
% No longer tries to figure out SPM version. MUST have SPM8 in path.
% Returns:
%   vol - SPM NIFTI volume structure (NOT the actual image data)
%   hdr - SPM Dicom file header structure
%------------------------------------------------------------------------

if (nargin < 1), infiles = {}; end
if (nargin < 2), outfile = ''; end
if (nargin < 3), make4D  = 0; end

default_outfile = '';

% --- prompt for dicom files ---
if (nargin < 1) || isempty(infiles)
    if (ispref(mfilename(),'dicompath')), dicompath = getpref(mfilename(),'dicompath');
    else, dicompath = [pwd() filesep()]; end
    [dcmfiles,dicompath] = uigetfile('*.dcm','Select the Dicom files to convert to NIFTI',dicompath,'MultiSelect','on');
    if (~iscell(dcmfiles) && ~ischar(dcmfiles)), return; end           % user hit cancel
    %dicompath = [dicompath filesep()];
    setpref(mfilename(),'dicompath',dicompath);                        % remember for next time
    infiles = strcat(dicompath,dcmfiles);
    pathpieces = strsplit(dicompath,filesep());
    default_outfile = pathpieces{end-1};                               % this get folder name holding dicoms
end

% --- prompt for NIFTI file name to save to ---
if (nargin < 2) || isempty(outfile)
    if (ispref(mfilename(),'niftipath')), niftipath = getpref(mfilename(),'niftipath');
    else, niftipath = [pwd() filesep()]; end
    [niftifile,niftipath] = uiputfile('*.nii','Enter the name of the output NIFTI file',[niftipath filesep() default_outfile]);
    if (~ischar(niftifile)), return; end                               % user hit cancel
    setpref(mfilename(),'niftipath',niftipath);                        % remember for next time
    outfile = [niftipath niftifile];
end

% --- Passed in cell array of filenames or dicominfo() structs? ---
if (iscell(infiles)), using_dicomhdrs = 0; else using_dicomhdrs = 1; end

% --- Build string array from cell array for SPM input ---
if (iscell(infiles))
    n = numel(infiles);
    for i=1:n
        if (using_dicomhdrs), tmp = infiles{i}.Filename;
        else                  tmp = infiles{i}; end
        
        if (i == 1), names = tmp;
        else         names = [names ; tmp]; end
    end
else
    if (using_dicomhdrs), names = infiles;
    else                  names = infiles.Filename; end
end

% --- Read dicom headers SPM-style ---
hdr = ME_spm8_dicom_headers3(names);

% --- Convert Dicoms to NIFTI using SPM8 functions ---
out = ME_spm8_dicom_convert3(hdr,'all',outfile);

% --- Read the Nifti volume structure(s) back from file(s) just created ---
vol = spm_vol(out.files);

% --- make a 4D NiIFTI if there were multiple volumes ---
if (numel(out.files) > 1)
    [outpath,outroot] = fileparts(outfile);
    filename4D = [outpath filesep() outroot '_4D.nii'];
    spm_file_merge(out.files,filename4D,0);  % SPM8
    matfile = [outpath filesep() outroot '_4D.mat'];
    delete(matfile);    % unwanted /mat file of xform matrices
end
return;	 
			 