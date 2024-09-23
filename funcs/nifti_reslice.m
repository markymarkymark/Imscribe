function [volout,outfile] = nifti_reslice(vol1,vol2,mat,newfilename,interp)

if nargin < 4, newfilename = '' ; end
if nargin < 5, interp = 4 ;       end  % interp = 0 means nearest neighbor (I THINK)

vols = [vol1 ; vol2]; % we're going to reslice the second volume to the first, using MAT
vols(2).mat = mat * vols(2).mat;

def_flags = struct('interp',interp,'mask',1,'mean',0,'which',1,'wrap',[0 0 0]','prefix','r');
outfile   = ME_spm8_reslice_images(vols,def_flags,newfilename);
volout    = spm_vol(outfile);
return

