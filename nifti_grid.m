% -------------------------------------------------------------------------
function [c,ctop,cbot,r,mat] = nifti_grid(volinfo)
% Returns coords of box around middle, top and bottom slice of a nifit volume
% Also returns (x,y,z) locations of all voxels in a Nifti Volume

use_spmstyle = 1;   % Nifti structure is from SPM read

% --- read SPM header style ---
if (use_spmstyle)
    nx  = volinfo.dim(1);
    ny  = volinfo.dim(2);
    nz  = volinfo.dim(3);
    x0 =  1;    x1 = nx+1;    
    y0 =  0;    y1 = ny;
    z0 =  0.5;  z1 = nz+0.5;
    mat = volinfo.mat(1:3,:);
    %mat = mat(1:3,:);       % remove 0, 0, 0, 1 row
    
% --- read load_untouch_nii() header style ---
else
    nx  = volinfo.hdr.dime.dim(2);
    ny  = volinfo.hdr.dime.dim(3);
    nz  = volinfo.hdr.dime.dim(4);
    x0 =  0;    x1 = nx-1;    
    y0 =  0;    y1 = ny-1;
    z0 =  0;    z1 = nz-1;
    mat = [volinfo.hdr.hist.srow_x ; volinfo.hdr.hist.srow_y ; volinfo.hdr.hist.srow_z];
end

% --- make corners of middle slice of volume ---
b4 = [x0 y0 (z1+z0)/2   1]';  
b3 = [x1 y0 (z1+z0)/2   1]';
b2 = [x1 y1 (z1+z0)/2   1]';
b1 = [x0 y1 (z1+z0)/2   1]';
b  = [b1' ; b2' ; b3' ; b4' ; b1'];	% add 1st point at end to close box end
c  = (mat * b')';  

% get box around top (of top) slice
b(:,3) = z1;
%b(:,3) = z1 + 0.5;
ctop   = (mat * b')';  

% get box around bottom (of bottom) slice
b(:,3) = z0;
%b(:,3) = z0 - 0.5;
cbot   = (mat * b')'; 

% --- generate x,y,z coords of all voxels in Nifti volume ---
if (nargout > 3)
    x0 = x0 + 0.5;  x1 = x1 - 0.5;
    y0 = y0 + 0.5;  y1 = y1 - 0.5;
    z0 = z0 + 0.5;  z1 = z1 - 0.5;
    
    nvox    = nx*ny*nz;
    [y,x,z] = meshgrid(y0:y1,x0:x1,z0:z1);       % Note X/Y switch needed fo meshgrid()!
    v       = zeros(4,nvox);                     % use 4 x nvox matrix for matmult later
    v(1,:)  = x(:);
    v(2,:)  = y(:);
    v(3,:)  = z(:);
    v(4,:)  = 1;
    clear x y z ;                                % free up memory
    r = mat * v;
    clear v ;
    r = reshape(r,3,nx,ny,nz);
end

return
