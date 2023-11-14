function [mat,rvol,outfile] = nifti_coreg(vol1,vol2,mode,newfilename,quality,quiet)
% Use SPM realign or affine coreg programs to coregisters TWO nifti volumes
%
% M.Elliott 2/07
% Check what version of Matlab to guess which SPM we are using
% mode = 'rigidbody' - intra-modal, rigid-body 6 param (i.e. SPM 'realign')
%      = 'affine'    - intra-modal, affine 12 param (i.e. SPM 'normalize')
%      = 'coreg'     - inter-modal, rigid-body 6 param (i.e. spm_coreg())
%
% NOTE: requires SPM8 in path!
%------------------------------------------------------------------------

if nargin < 4, newfilename = '' ; end
if nargin < 5, quality     = 1  ; end
if nargin < 5, quiet       = 0  ; end

mat     = [];
rvol    = [];
outfile = '';

% --- Do rigid-body (i.e. intra-subject) "realignment" ---
switch (mode)
    case 'rigidbody_fast'
        flags = struct('quality',quality,'fwhm',5,'sep',4,'interp',2,'wrap',[0 0 0],'rtm',0,'PW','','graphics',0,'lkp',1:6);
        
        % --- only the reference volume provided, initialize for faster performance ---
        if (isempty(vol2))
        	ME_realign_series_fast(vol1,flags);
            return
        end
        
        % --- Do rigid body realign, using initialized reference volume from above --
        vols  = [vol1 ; vol2];	% Make array of vol structures and call SPM2's realign_series()
        vols2 = ME_realign_series_fast(vols,flags);
        
    case 'rigidbody'
        if (~quiet), disp('Performing rigid-body INTRA-modal coregistration...'); end
        vols  = [vol1 ; vol2];	% Make array of vol structures and call SPM2's realign_series()
        flags = struct('quality',0.75,'fwhm',5,'sep',4,'interp',2,'wrap',[0 0 0],'rtm',0,'PW','','graphics',0,'lkp',1:6);
        vols2 = ME_realign_series(vols,flags);
        
    case 'affine'
        if (~quiet), disp('Performing affine INTRA-modal coregistration...'); end
        flags = struct('smosrc',6,'smoref',6,'regtype','subj','cutoff',30,'nits',16,'reg',0.1,'graphics',0);
        pnorm = ME_spm8_normalise(vol1,vol2,'','','',flags);
        vols2 = [vol1 ; vol2];
        vols2(2).mat = pnorm.M * vols2(2).mat;		% incorporate xform into target volume for reslice below
        
    case {'coreg','intermodal','intramodal'}  % note bug in 'intramodal' option!
        if (~quiet), disp('Performing rigid-body INTER-modal coregistration...'); end
%        pmat  = ME_spm_coreg(vol1,vol2);
        pmat  = spm_coreg(vol1,vol2);
        M     = inv(spm_matrix(pmat));
        vols2 = [vol1 ; vol2];
        vols2(2).mat = M * vols2(2).mat;		% incorporate xform into target volume for reslice below
        
    otherwise
        error('Unknown mode option\n');
end

%  --- Return transform matrix that maps TEMPLATE to TARGET ---
%mat = vols2(1).mat/vols2(2).mat;
%mat = inv(pnorm.M);
mat = inv(vols2(2).mat/vol2.mat); % This is the same as pnorm.M, but works for regid-body as well.

% --- Compute coreg'd image (i.e. SPM calls it "resliced") ---_
if (isequal(mode,'rigidbody_fast'))
    if (isempty(newfilename)), return; end  % don't need resliced image if only using it for RTfmri monitoring of motion
    flags = struct('interp',4,'mask',0,'mean',0,'which',1,'wrap',[0 0 0]','prefix','r');
else
    flags = struct('interp',4,'mask',1,'mean',0,'which',1,'wrap',[0 0 0]','prefix','r');
end    

% --- Return the regridded & coreg'd NIFTI header struct ---
outfile = ME_spm8_reslice_images(vols2,flags,newfilename);
rvol = spm_vol(outfile);
return;	 
			 