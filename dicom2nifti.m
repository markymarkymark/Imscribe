function [vol,hdr,filename4D] = dicom2nifti(infiles,outfile)
% Use SPM to convert Dicoms to Nifti volume.
% Requires SPM8 to be in the MATLAB path
%
% Syntax:
%   [vol,hdr,filename4D] = dicom2nifti([infiles],[outfile])
% Inputs:
%   infiles - a cell array of dicom filenames, OR cell array of dicominfo() structs
%           - if missing or empty, user is prompted to choose the files
%   outfile - filename root for result (*.nii)
%           - if missing or empty, user is prompted to specify
% Outputs:
%   vol     - SPM NIFTI volume structure (NOT the actual image data)
%   hdr     - SPM Dicom file header structure
%   filename4D - filename for the created 4D Nifti (if input dicoms are multiple volumes)
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

vol = [];
if (nargin < 1), infiles = {}; end
if (nargin < 2), outfile = ''; end

SPMver = get_spm_version('SPM8');  % MELLIOTT
if (isempty(SPMver)), return; end

% --- prompt for dicom files ---
default_outfile = '';
if (isempty(infiles))
    [infiles,~,~,~,~,botdir] = uigetfile_plus(mfilename(),'dicompath',pwd(),{'*.dcm;*.MR;*.IMA';'*.*'},'Select the Dicom files to convert to NIFTI','MultiSelect','on');
    if (isempty(infiles)), return; end                                   % user hit cancel
    default_outfile = botdir{1}(1:end-1);                                % this gets folder name holding dicoms, removing final '/'
end

% --- prompt for NIFTI file name to save to ---
if (isempty(outfile))
    outfile = uiputfile_plus(mfilename(),'niftipath',pwd(),'*.nii','Enter the name of the output NIFTI file',default_outfile);
    if (isempty(outfile)), return; end                                   % user hit cancel
end

% --- Passed in cell array of filenames or dicominfo() structs? ---
if (iscell(infiles)), using_dicomhdrs = 0; else, using_dicomhdrs = 1; end

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
hdr = me_spm8_dicom_headers(names);

% --- Convert Dicoms to NIFTI using SPM8 functions ---
out = me_spm8_dicom_convert(hdr,'all',outfile);

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
			 