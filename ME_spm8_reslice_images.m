%_______________________________________________________________________
function [outfilename] = ME_spm8_reslice_images(P,flags,newfilename)
% Reslices images volume by volume
% FORMAT reslice_images(P,flags)
%
% P        - matrix of image handles from spm_vol.
%            All operations are performed relative to the first image.
%            ie. resampling of images is into the space of the first image.
%
% flags    - a structure containing various options.  The fields are:
%
%         mask - mask output images (1 for yes, 0 for no)
%                To avoid artifactual movement-related variance the realigned
%                set of images can be internally masked, within the set (i.e.
%                if any image has a zero value at a voxel than all images have
%                zero values at that voxel).  Zero values occur when regions
%                'outside' the image are moved 'inside' the image during
%                realignment.
%
%         mean - write mean image
%                The average of all the realigned scans is written to
%                mean*.img.
%
%         interp - the B-spline interpolation method (see spm_bsplinc and spm_bsplins).
%                Non-finite values result in Fourier interpolation
%
%         which - Values of 0, 1 or 2 are allowed.
%                0   - don't create any resliced images.
%                      Useful if you only want a mean resliced image.
%                1   - don't reslice the first image.
%                      The first image is not actually moved, so it may not be
%                      necessary to resample it.
%                2   - reslice all the images.
%         wrap - three values of either 0 or 1, representing wrapping in each of
%                the dimensions.  For fMRI, [1 1 0] would be used.  For PET, it would
%                be [0 0 0].
%
%         prefix - prefix for resliced images. Defaults to 'r'.
%
%             The spatially realigned images are written to the orginal
%             subdirectory with the same filename but prefixed with an 'r'.
%             They are all aligned with the first.
%
% Modified by MElliott 
%   - cut from spm_reslice.m
%   - User provides output filename. 
%       Hence, only works for one output volume.
%       Hence, should always be called with flags.which = 1
%   - forced removal of NaNs from result



% MElliott
if nargin < 3, newfilename = '' ; end

if ~isfinite(flags.interp), % Use Fourier method
    % Check for non-rigid transformations in the matrixes
    for i=1:prod(size(P)),
        pp = P(1).mat\P(i).mat;
        if any(abs(svd(pp(1:3,1:3))-1)>1e-7),
            fprintf('\n  Zooms  or shears  appear to  be needed');
            fprintf('\n  (probably due to non-isotropic voxels).');
            fprintf('\n  These  can not yet be  done  using  the');
            fprintf('\n  Fourier reslicing method.  Switching to');
            fprintf('\n  7th degree B-spline interpolation instead.\n\n');
            flags.interp = 7;
            break;
        end;
    end;
end;

if flags.mask | flags.mean,
    x1    = repmat((1:P(1).dim(1))',1,P(1).dim(2));
    x2    = repmat( 1:P(1).dim(2)  ,P(1).dim(1),1);
    if flags.mean,
        Count    = zeros(P(1).dim(1:3));
        Integral = zeros(P(1).dim(1:3));
    end;
    if flags.mask, msk = cell(P(1).dim(3),1);  end;
    for x3 = 1:P(1).dim(3),
        tmp = zeros(P(1).dim(1:2));
        for i = 1:prod(size(P)),
            tmp = tmp + getmask(inv(P(1).mat\P(i).mat),x1,x2,x3,P(i).dim(1:3),flags.wrap);
        end;
        if flags.mask, msk{x3} = find(tmp ~= prod(size(P))); end;
        if flags.mean, Count(:,:,x3) = tmp; end;
    end;
end;

nread = prod(size(P));
if ~flags.mean,
    if flags.which == 1, nread = nread - 1; end;
    if flags.which == 0, nread = 0; end;
end;

tiny = 5e-2; % From spm_vol_utils.c

[x1,x2] = ndgrid(1:P(1).dim(1),1:P(1).dim(2));

PO    = P;
nread = 0;
d     = [flags.interp*[1 1 1]' flags.wrap(:)];

for i = 1:prod(size(P)),

    if (i>1 & flags.which==1) | flags.which==2, write_vol = 1; else, write_vol = 0; end;
    if write_vol | flags.mean,                   read_vol = 1; else   read_vol = 0; end;

    if read_vol,
        if ~isfinite(flags.interp),
            v = abs(kspace3d(spm_bsplinc(P(i),[0 0 0 ; 0 0 0]'),P(1).mat\P(i).mat));
            for x3 = 1:P(1).dim(3),
                if flags.mean,
                    Integral(:,:,x3) = Integral(:,:,x3) + ...
                        nan2zero(v(:,:,x3).*getmask(inv(P(1).mat\P(i).mat),x1,x2,x3,P(i).dim(1:3),flags.wrap));
                end;
                if flags.mask, tmp = v(:,:,x3); tmp(msk{x3}) = NaN; v(:,:,x3) = tmp; end;
            end;
        else,
            C = spm_bsplinc(P(i), d);
            v = zeros(P(1).dim);
            for x3 = 1:P(1).dim(3),

                [tmp,y1,y2,y3] = getmask(inv(P(1).mat\P(i).mat),x1,x2,x3,P(i).dim(1:3),flags.wrap);
                v(:,:,x3)      = spm_bsplins(C, y1,y2,y3, d);
                % v(~tmp)        = 0;

                if flags.mean, Integral(:,:,x3) = Integral(:,:,x3) + nan2zero(v(:,:,x3)); end;

                if flags.mask, tmp = v(:,:,x3); tmp(msk{x3}) = NaN; v(:,:,x3) = tmp; end;
            end;
        end;
        if write_vol,
            VO         = P(i);
%            [pth,nm,xt,vr] = fileparts(deblank(P(i).fname));
%            VO.fname   = fullfile(pth,[flags.prefix nm xt vr]);
            if (isempty(newfilename))                               % MELLIOTT
                [pth,nm,xt] = fileparts(deblank(P(i).fname));       % Removed 'vr' MELLIOTT
                outfilename = fullfile(pth,[flags.prefix nm xt]);
            else
                outfilename = fullfile(newfilename);
            end
            VO.fname    = outfilename;    
            VO.dim     = P(1).dim(1:3);
            VO.dt      = P(i).dt;
            VO.pinfo   = P(i).pinfo;
            VO.mat     = P(1).mat;
            VO.descrip = 'spm - realigned';
v(isnan(v)) = 0;    % MELLIOTT
            VO = spm_write_vol(VO,v);
        end;

        nread = nread + 1;
    end;
end;

if flags.mean
    % Write integral image (16 bit signed)
    %-----------------------------------------------------------
    Integral   = Integral./Count;
    PO         = P(1);
    PO         = rmfield(PO,'pinfo');
    [pth,nm,xt,vr] = fileparts(deblank(P(1).fname));
    PO.fname   = fullfile(pth,['mean' nm xt]);
    PO.pinfo   = [max(max(max(Integral)))/32767 0 0]';
    PO.descrip = 'spm - mean image';
    PO.dt      = [4 spm_platform('bigend')];
    spm_write_vol(PO,Integral);
end

return;
%_______________________________________________________________________

%_______________________________________________________________________
function v = kspace3d(v,M)
% 3D rigid body transformation performed as shears in 1D Fourier space.
% FORMAT v1 = kspace3d(v,M)
% Inputs:
% v - the image stored as a 3D array.
% M - the rigid body transformation matrix.
% Output:
% v - the transformed image.
%
% The routine is based on the excellent papers:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI
% Submitted to MRM (April 1999) and avaliable from:
% http://varda.biophysics.mcw.edu/~cox/index.html.
% and:
% W. F. Eddy, M. Fitzgerald and D. C. Noll (1996)
% Improved Image Registration by Using Fourier Interpolation
% Magnetic Resonance in Medicine 36(6):923-931
%_______________________________________________________________________

[S0,S1,S2,S3] = shear_decomp(M);

d  = [size(v) 1 1 1];
g = 2.^ceil(log2(d));
if any(g~=d),
    tmp = v;
    v   = zeros(g);
    v(1:d(1),1:d(2),1:d(3)) = tmp;
    clear tmp;
end;

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2),
    t        = reshape( exp((j*S3(3,2) + S3(3,1)*(1:g(1)) + S3(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
    v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end;

% XZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(2)-1)/2) 0 (-g(2)/2+1):-1])/g(2);
for k=1:g(3),
    t        = exp( (k*S2(2,3) + S2(2,1)*(1:g(1)) + S2(2,4)).'*tmp1);
    v(:,:,k) = real(ifft(fft(v(:,:,k),[],2).*t,[],2));
end;

% YZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(1)-1)/2) 0 (-g(1)/2+1):-1])/g(1);
for k=1:g(3),
    t        = exp( tmp1.'*(k*S1(1,3) + S1(1,2)*(1:g(2)) + S1(1,4)));
    v(:,:,k) = real(ifft(fft(v(:,:,k),[],1).*t,[],1));
end;

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2),
    t        = reshape( exp( (j*S0(3,2) + S0(3,1)*(1:g(1)) + S0(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
    v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end;

if any(g~=d), v = v(1:d(1),1:d(2),1:d(3)); end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [S0,S1,S2,S3] = shear_decomp(A)
% Decompose rotation and translation matrix A into shears S0, S1, S2 and
% S3, such that A = S0*S1*S2*S3.  The original procedure is documented
% in:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI

A0 = A(1:3,1:3);
if any(abs(svd(A0)-1)>1e-7), error('Can''t decompose matrix'); end;


t  = A0(2,3); if t==0, t=eps; end;
a0 = pinv(A0([1 2],[2 3])')*[(A0(3,2)-(A0(2,2)-1)/t) (A0(3,3)-1)]';
S0 = [1 0 0; 0 1 0; a0(1) a0(2) 1];
A1 = S0\A0;  a1 = pinv(A1([2 3],[2 3])')*A1(1,[2 3])';  S1 = [1 a1(1) a1(2); 0 1 0; 0 0 1];
A2 = S1\A1;  a2 = pinv(A2([1 3],[1 3])')*A2(2,[1 3])';  S2 = [1 0 0; a2(1) 1 a2(2); 0 0 1];
A3 = S2\A2;  a3 = pinv(A3([1 2],[1 2])')*A3(3,[1 2])';  S3 = [1 0 0; 0 1 0; a3(1) a3(2) 1];

s3 = A(3,4)-a0(1)*A(1,4)-a0(2)*A(2,4);
s1 = A(1,4)-a1(1)*A(2,4);
s2 = A(2,4);
S0 = [[S0 [0  0 s3]'];[0 0 0 1]];
S1 = [[S1 [s1 0  0]'];[0 0 0 1]];
S2 = [[S2 [0 s2  0]'];[0 0 0 1]];
S3 = [[S3 [0  0  0]'];[0 0 0 1]];
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
tiny = 5e-2; % From spm_vol_utils.c
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
Mask = logical(ones(size(y1)));
if ~wrp(1), Mask = Mask & (y1 >= (1-tiny) & y1 <= (dim(1)+tiny)); end;
if ~wrp(2), Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny)); end;
if ~wrp(3), Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny)); end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vo = nan2zero(vi)
vo = vi;
vo(~isfinite(vo)) = 0;
return;
%_______________________________________________________________________

%_______________________________________________________________________

