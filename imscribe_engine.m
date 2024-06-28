function [vols,cboxes,mat,opt] = imscribe_engine(files1,files2,files3,opt)

% --- Set return variables ---
vols = {}; cboxes = {}; mat = []; 
if (nargin < 4), return; end
opt.useXform       = checkstruct(opt,'useXform',[]);   
opt.winid          = checkstruct(opt,'winid',[]);   
opt.txtid          = checkstruct(opt,'txtid',[]);   
opt.scratchpath    = checkstruct(opt,'scratchpath',pwd());   
opt.svsmask        = checkstruct(opt,'svsmask',0);   
opt.drawsourceonly = checkstruct(opt,'drawsourceonly',0);   
opt.contourlevel   = checkstruct(opt,'contourlevel',50);   
opt.slicenum       = checkstruct(opt,'slicenum',50);   
opt.imbright       = checkstruct(opt,'imbright',1);   
opt.displane       = checkstruct(opt,'displane',1);   
opt.coregmethod    = checkstruct(opt,'coregmethod','rigidbody');   
opt.datatypes      = checkstruct(opt,'datatypes',{'NIFTI','NIFTI','NIFTI'});   % assume NIFTIs

% --- Set some local variables ---
mat_N2S        = [-1 0 0; 0 -1 0 ; 0 0 1]'; % xform from NIFTI to Siemens MRI coord space ---
mat_S2N        = inv(mat_N2S);

tic;
disp(' ');
disp('---------------------------------------------------------------------------------');
fprintf(1,'Starting %s...\n',mfilename());
fprintf(1,'\tSOURCE Structural: %s\n',files1{1});
fprintf(1,'\tSOURCE FOV/ROI   : %s\n',files2{1});
if (~opt.svsmask && ~opt.drawsourceonly)
    fprintf(1,'\tTARGET Structural: %s\n',files3{1});
end

disp('---------------------------------------------------------------------------------');
disp('Reading data files...');
%--- SOURCE STRUCT ---
niftifile1 = [opt.scratchpath 'SourceStruct.nii'];
switch opt.datatypes{1}
    case 'DICOM'
        fprintf(1,'\tConverting %1d SOURCE Structural Dicoms to Nifti\n',numel(files1));
        vol1 = dicom2nifti(files1,niftifile1);
        if (iscell(vol1)), vol1 = vol1{1}; end
    case 'NIFTI'
        fprintf(1,'\tReading SOURCE Structural Nifti header\n');
        vol1 = spm_vol(files1{1});
        copyfile(files1{1},niftifile1,'f');
        fprintf(1,'\tSaved a copy of the SOURCE Structural to %s\n',niftifile1);
    otherwise
        fprintf(2,'ERROR: Unrecognized data file type "%s"\n',opt.datatypes{1}); return; 
end
c1    = nifti_grid(vol1,1,1);     % get bounding box of vol1 FOV
c1    = box_prescrip(c1,mat_N2S); % get Siemens prescription of vol1
I1vol = spm_read_vols(vol1);

% --- SOURCE FOV/SVS ---
niftifile2 = [opt.scratchpath 'SourceFOV.nii'];
switch opt.datatypes{2}
    case 'RDA'
        is_spec   = 1;
        fprintf(1,'\tCalculating RDA SVS location\n');
        [~,rdahdr] = readRDA(files2{1});
        c2         = rda_3dcalc(rdahdr,mat_S2N);
    case 'DICOM'
        dcmhdr    = dicominfo(files2{1},'UseDictionaryVR',true);
        is_spec   = isfield(dcmhdr,'Private_7fe1_1010') || isfield(dcmhdr,'SpectroscopyData');
        is_mosaic = ~isempty(strfind(dcmhdr.ImageType,'MOSAIC'));
        if ((is_spec || is_mosaic) && numel(files2) > 1)
            fprintf(2,'\tWARNING: Multiple Spectro or Mosaic Dicom files selected for SOURCE SVS/FOV. Using first file only.\n');
            files2 = {files2{1}};
        end
        if (~is_spec)
            fprintf(1,'\tConverting %1d SOURCE FOV Dicoms to Nifti\n',numel(files2));
            vol2 = dicom2nifti(files2,niftifile2);
            if (iscell(vol2)), vol2 = vol2{1}; end
            c2 = nifti_grid(vol2,0,1); 
        else
            fprintf(1,'\tCalculating Dicom SVS location\n');
            c2 = dicom_3dcalc(dcmhdr,[],[],mat_S2N);
            is_roi = 0; % ignore ROI threshold for SVS
        end
    case 'NIFTI'
        is_spec = 0;
        fprintf(1,'\tReading SOURCE FOV/ROI Nifti header\n');
        vol2 = spm_vol(files2{1});
        if (numel(vol2) > 1), vol2 = vol2(1); end   % in case 4D Nifti was chosen
        c2 = nifti_grid(vol2,0,1); 
        copyfile(files2{1},niftifile2,'f');
        fprintf(1,'\tSaved a copy of the SOURCE FOV/ROI to %s\n',niftifile2);
    otherwise
        fprintf(2,'ERROR: Unrecognized data file type "%s"\n',opt.datatypes{2}); return; 
end
c2 = box_prescrip(c2,mat_N2S);
%I2vol = spm_read_vols(vol2);  
I2vol = [];   % don't actually need voxel data for FOV/SVS

%--- TARGET STRUCT ---
if (~opt.svsmask && ~opt.drawsourceonly)
    niftifile3 = [opt.scratchpath 'TargetStruct.nii'];
    switch opt.datatypes{3}
        case 'DICOM'
            fprintf(1,'\tConverting %1d TARGET Structural Dicoms to Nifti\n',numel(files3));
            vol3 = dicom2nifti(files3,niftifile3);
            if (iscell(vol3)), vol3 = vol3{1}; end
        case 'NIFTI'
            fprintf(1,'\tReading TARGET Structural Nifti header\n');
            vol3 = spm_vol(files3{1});
            copyfile(files3{1},niftifile3,'f');
            fprintf(1,'\tSaved a copy of the TARGET Structural to %s\n',niftifile3);
        otherwise
            fprintf(2,'ERROR: Unrecognized data file type "%s"\n',opt.datatypes{3}); return; 
    end
    c3    = nifti_grid(vol3,1,1); 
    c3    = box_prescrip(c3,mat_N2S);
    I3vol = spm_read_vols(vol3); 

else
    c3    = [];
    I3vol = [];
end

disp('---------------------------------------------------------------------------------');
disp('Siemens-style Prescriptions...');
fprintf(1,'\tSOURCE Struct: %s \n',c1.prescrip);
fprintf(1,'\tTEMPLATE FOV:  %s (extent = %1.1f x %1.1f x %1.1f)\n',c2.prescrip,c2.lx,c2.ly,c2.lz);
if (~opt.svsmask && ~opt.drawsourceonly), fprintf(1,'\tTARGET Struct: %s\n',c3.prescrip); end

% --- Display initial config ---
if (~isempty(opt.winid))
    imscribe_draw(I1vol,I3vol,I1vol,c1,c3,c1,c2,opt,[2 1 1]);
end

% --- Only draw Template LOC/FOV overlay, thenn exit ---
if (opt.drawsourceonly)
    fprintf(1,'Only asked to draw overlay of Source Struct & FOV. Done.\n');
    vols   = {I1vol, I2vol, I3vol, []};
    cboxes = {c1, c2, c3, [], []};
    return
end

% --- Generate NIFTI with SVS or FOV mask on SOURCE Struct grid ---
if (opt.svsmask)
    if (is_spec)
        disp('---------------------------------------------------------------------------------');
        disp('Generating SVS mask in Template Structural space...');
        fname = [opt.scratchpath 'svsmask.nii'];
    else
        disp('---------------------------------------------------------------------------------');
        disp('Generating FOV mask in Template Structural space...');
        fname = [opt.scratchpath 'fovmask.nii'];
    end
    if (~contains(c2.prescrip,'<') || (~contains(c2.prescrip,'>') || abs(c2.rot) > 1))
        fprintf(1,'Handling oblique FOV/SVS voxel...\n');
        [v1,v2,v3] = box_orient(c2.cmid);
        Vsvs_in =  [v1 v2 v3];
        Vsvs_out = [[1;0;0] [0;1;0] [0;0;1]];
        Msvs = Vsvs_out/Vsvs_in;     % this is the transform that rotates the SVS voxel into a non-oblique orientation
        cbot2x = (Msvs * c2.cbot')';
        ctop2x = (Msvs * c2.ctop')';
        r1volx =  Msvs * c1.rvol(:,:);  % now the voxel coords of the Template Loc are rotated so the SVS voxel is non-oblique
    else
        ctop2x = c2.ctop;
        cbot2x = c2.cbot;
        r1volx = c2.rvol;
    end
    
    % make 3D image with mask of pixels inside SVS FOV
    Isvs = uint16(I1vol * 0);
    xr = [min([cbot2x(:,1) ;ctop2x(:,1)]) max([cbot2x(:,1) ;ctop2x(:,1)])];
    yr = [min([cbot2x(:,2) ;ctop2x(:,2)]) max([cbot2x(:,2) ;ctop2x(:,2)])];
    zr = [min([cbot2x(:,3) ;ctop2x(:,3)]) max([cbot2x(:,3) ;ctop2x(:,3)])];
    psvs =  (r1volx(1,:) >= xr(1)) & (r1volx(1,:) <= xr(2)) & ...
            (r1volx(2,:) >= yr(1)) & (r1volx(2,:) <= yr(2)) & ...
            (r1volx(3,:) >= zr(1)) & (r1volx(3,:) <= zr(2)) ;
    if (~any(psvs))
        fprintf(2,'ERROR: No voxels in SOURCE Structural fall inside the FOV/SVS ROI! No mask image generated.\n');
        return
    end
    Isvs(psvs) = 1;
    
    % write the Nifti file with SVS mask
    svsvol = vol1;   % copy SPM image volume of SOURCE Struct
    svsvol.fname = fname;
    svsvol.descrip = 'FOV/SVS mask';
    spm_write_vol(svsvol,Isvs);
    fprintf(1,'Saved FOV/SVS mask in %s\n',svsvol.fname);
    
    vols   = {I1vol, I2vol, I3vol, []};
    cboxes = {c1, c2, c3, [], []};
    return
end

% --- Register SOURCE Struct to TARGET Struct ---
niftifile4 = [opt.scratchpath 'Target2Source.nii'];
if (isempty(opt.useXform))
    disp('---------------------------------------------------------------------------------');
    disp('Registering TARGET Struct to SOURCE Struct...');
    [mat,vol4] = nifti_coreg(vol1,vol3,opt.coregmethod,niftifile4); % this returns re-sliced vol3 on vol1 grid
    fprintf(1,'Saved result in %s\n',vol4.fname);
else
    disp('---------------------------------------------------------------------------------');
    disp('Using provided xform betwen TARGET and SOURCE Structurals...');
    mat = opt.useXform;
%    vol4 = nifti_reslice(vol1,vol3,mat,niftifile4);
    vol4 = nifti_reslice(vol3,vol1,mat,niftifile4);
end
c4    = nifti_grid(vol4,1,1);     % get bounding box of vol1 FOV
c4    = box_prescrip(c4,mat_N2S); % get Siemens prescription of vol1
I4vol = spm_read_vols(vol4);

% --- Reslice the SOURCE FOV into TARGET Struct space (for off-line checking ONLY - not used here) ---
% if (~is_spec)
%     disp('---------------------------------------------------------------------------------');
%     disp('Transforming SOURCE FOV to TARGET Struct space...');
%     niftifile5 = [opt.scratchpath 'SourceFOV2Target.nii'];
%     rvol2 = nifti_reslice(vol3,vol2,mat,niftifile5);
%     fprintf(1,'Saved result in %s\n',rvol2.fname);
% end

% --- Apply transform to bounding box of FOV/SVS ---
c2x      = c2;
c2x.cmid = transform_coords(c2x.cmid,mat);
c2x.ctop = transform_coords(c2x.ctop,mat);
c2x.cbot = transform_coords(c2x.cbot,mat);
% c2x.cmid = transform_coords(c2x.cmid,inv(mat));
% c2x.ctop = transform_coords(c2x.ctop,inv(mat));
% c2x.cbot = transform_coords(c2x.cbot,inv(mat));
c2x      = box_prescrip(c2x,mat_N2S);

imscribe_draw(I1vol,I4vol,I4vol,c1,c4,c4,c2x,opt,[4 3 2]);

% --- Display the coreg matrix params ---
disp('---------------------------------------------------------------------------------');
disp('Coregistration params...');
params  = spm_imatrix(mat);
fprintf(1,'Tx,Ty,Tz = %6.2f,%6.2f,%6.2f\n',params(1),params(2),params(3));
fprintf(1,'R1,R2,R3 = %6.2f,%6.2f,%6.2f\n',rad2deg(params(4)),rad2deg(params(5)),rad2deg(params(6)));
fprintf(1,'Sx,Sy,Sz = %6.2f,%6.2f,%6.2f\n',params(7),params(8),params(9));
fprintf(1,'Ax,Ay,Az = %6.2f,%6.2f,%6.2f\n',params(10),params(11),params(12));
disp('---------------------------------------------------------------------------------');
fprintf(1,'New FOV prescription: %s (rot = %1.1f)\n',c2x.prescrip,c2x.rot);
fprintf(1,'\t           extent: %1.1f x %1.1f x %1.1f\n',c2x.lx,c2x.ly,c2x.lz);
disp('---------------------------------------------------------------------------------');

% --- Write new prescription to a text file, appending if exists ---
pfile = [opt.scratchpath filesep() 'imscribe_prescriptions.txt'];
fp=fopen(pfile,'a+');
if (fp ~= -1)
    fprintf(fp,'%s   %s (rot = %1.1f)\n',datestr(now()),c2x.prescrip,c2x.rot);
    fclose(fp);
end
fprintf(1,'Wrote prescription to %s\n',pfile)

disp('---------------------------------------------------------------------------------');
disp('---------------------------------------------------------------------------------');
disp('Done.');
toc;

vols   = {I1vol, I2vol, I3vol, I4vol};
cboxes = {c1, c2, c3, c4, c2x};
end


% -------------------------------------------------------------------------
% SUBROUTINES
% -------------------------------------------------------------------------

%------------------------------------------------------------------------
function res = box_prescrip(res,xform)
% Return string with Siemens-style prescription
%	e.q. L10 P5 H12.5 T>C-20
% Input is bounding box coords, as returned by dicom_3dcalc(), or nifti_grid()
%   NOTE: the bounding box should be in Siemens coord space (i.e. nifti_grid() result needs mat_N2S*)
%------------------------------------------------------------------------
cmid = res.cmid;
ctop = res.ctop;
cbot = res.cbot;
if (~isempty(xform))
    cmid = (xform*cmid')';
    ctop = (xform*ctop')';
    cbot = (xform*cbot')';
end

% --- get inplane & normal vectors, and origin of bounding box ---
[v1,v2,v3,res.origin] = box_orient(cmid);

% --- Get Siemens desciption of slice plane orientation ---
res.plane = vector_prescrip(v3);

% --- Add center voxel location ---
res.center = origin_prescrip(res.origin);
res.prescrip = [res.plane ' ' res.center];

% --- Determine the inplane rotation ---
res.rot = inplane_rot(res.prescrip,v1,v2,v3);

% extent
res.lx = sqrt( (cmid(2,1)-cmid(1,1))^2 + (cmid(2,2)-cmid(1,2))^2 + (cmid(2,3)-cmid(1,3))^2 );
res.ly = sqrt( (cmid(3,1)-cmid(2,1))^2 + (cmid(3,2)-cmid(2,2))^2 + (cmid(3,3)-cmid(2,3))^2 );
res.lz = sqrt( (ctop(1,1)-cbot(1,1))^2 + (ctop(1,2)-cbot(1,2))^2 + (ctop(1,3)-cbot(1,3))^2 );
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
function cout = transform_coords(cin,mat)
% xform XYZ coords
ctmp = [cin ones(5,1)];     % make each XYZ into XYZ1
cout = mat * ctmp';			% transform 
cout = cout(1:3,:)';		% clip XYZ1 into XYZ
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function rot = inplane_rot(prescrip,v1,v2,v3)
% figure out in-plane rotation

X = [1 0 0]'; % cartesian unit vectors
Y = [0 1 0]'; 
Z = [0 0 1]'; 

% get plane and build rotation matrix that produces this plane's orientation
plane = prescrip(1:3);
p     = strfind(prescrip,'>');
nangles = numel(p);
switch(nangles)
    case 0
        angle1 = 0;
        angle2 = 0;
    case 1
        num = sscanf(prescrip,'%*c>%*c%f %*s');
        angle1 = num;
        angle2 = 0;
    case 2
        nums = sscanf(prescrip,'%*c>%*c%f>%c%f %*s');
        angle1 = nums(1);
        plane  = [plane '>' char(nums(2))]; % added 2nd tilt plane
        angle2 = nums(3);
end

rotsign = -1;
switch plane
	case 'TRA'
        R = rotx(0.);               % Transverse slice is unrotated coord system
	case 'COR'
        R = rotx(deg2rad(-90),1);    
	case 'SAG'
        R = rotx(deg2rad(-90),1) * rotx(deg2rad(-90),2);   
        
    case 'T>C'
        R = rotx(deg2rad(angle1),1);   
	case 'T>S'
        R = rotx(deg2rad(-angle1),2);  
        
    case 'C>S'
        R1 = rotx(deg2rad(-90),1);    
        R2 = rotx(deg2rad(angle1),3);  
        R  = R2*R1;
    case 'C>T'
        R1 = rotx(deg2rad(-90),1);    
        R2 = rotx(deg2rad(-angle1),1);  
        R  = R2*R1;

    case 'S>C'
        R1 = rotx(deg2rad(-90),1) * rotx(deg2rad(-90),2);    
        R2 = rotx(deg2rad(-angle1),3);
        R  = R2*R1;
    case 'S>T'
        R1 = rotx(deg2rad(-90),1) * rotx(deg2rad(-90),2);    
        R2 = rotx(deg2rad(angle1),2);
        R  = R2*R1;
        rotsign = 1;
        
    case 'T>S>C'
        R1 = rotx(deg2rad(-angle1),2); 
        R2 = rotx(deg2rad(angle2),1); 
        R  = R2*R1;
    case 'T>C>S'    % This may not be right??
        R1 = rotx(deg2rad(angle1),1);   
        R2 = rotx(deg2rad(-angle2),2); 
        R  = R2*R1;

    case 'C>T>S'
        R1 = rotx(deg2rad(-90),1);    
        R2 = rotx(deg2rad(-angle1),1);  
        R3 = rotx(deg2rad(angle2),3);  
        R  = R3*R2*R1;
    case 'C>S>T'
        R1 = rotx(deg2rad(-90),1);    
        R2 = rotx(deg2rad(angle1),3);  
        R3 = rotx(deg2rad(angle2),1);  
        R  = R3*R2*R1;
        
    case 'S>T>C'
        R1 = rotx(deg2rad(-90),1) * rotx(deg2rad(-90),2);    
        R2 = rotx(deg2rad(angle1),2);
        R3 = rotx(deg2rad(-angle2),3);
        R  = R3*R2*R1;
        rotsign = 1;
    case 'S>C>T'  % This may not be right??
        R1 = rotx(deg2rad(-90),1) * rotx(deg2rad(-90),2);    
        R2 = rotx(deg2rad(-angle1),3);
        R3 = rotx(deg2rad(angle2),2);
        R  = R3*R2*R1;
              
end       

% get expected row and column direction vectors if there were no inplane rotation      
v1x = R*X;  
v2x = R*Y;  

% undo prescrip rotation to get v1,v2 vectors in Transverse slice
Rinv = inv(R);
v1p = Rinv*v1;  % should be [1,0,0] if no inplane rotation
v2p = Rinv*v2;  % should be [0,1,0] ...
rot = rotsign * rad2deg(atan2(v1p(2),v1p(1)));
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
function trans = rotx(theta,axis,use_4by4)

if (nargin < 2), axis = 1; end % default is rotate about X
if (nargin < 3), use_4by4 = 0; end % default is generate only 3x3 Euler

% compute coefs
cot = cos(theta);
sit = sin(theta);

% create transformation
switch axis
    case 1
    trans = [...
        1 0 0 0;...
        0 cot -sit 0;...
        0 sit cot 0;...
        0 0 0 1];
    
    case 2    
    trans = [...
        cot  0  sit  0;...
        0    1    0  0;...
        -sit 0  cot  0;...
        0    0    0  1];
    
    case 3
    trans = [...
        cot -sit 0 0;...
        sit  cot 0 0;...
        0 0 1 0;...
        0 0 0 1];
end

if (~use_4by4)
    trans = trans(1:3,1:3);
end
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function [v1,v2,v3,origin] = box_orient(c)
% calc description of fov box given by it's corners.
%	returns:
%       v1,v2,v3:   row,column,normal vectors  
%       origin:     XYZ coord of center
%
% M.Elliott 2/07
%------------------------------------------------------------------------

% --- calc V1, V2, V3 - the row, column, and normal vectors of the image plane
v1 = [c(2,1)-c(1,1) ; c(2,2)-c(1,2) ; c(2,3)-c(1,3) ];	% increasing columns
v2 = [c(3,1)-c(2,1) ; c(3,2)-c(2,2) ; c(3,3)-c(2,3) ];	% increasing rows
v3 = cross(v1,v2);

% --- normalize directions vectors ---
v1 = v1/sqrt(sum(v1.^2));
v2 = v2/sqrt(sum(v2.^2));
v3 = v3/sqrt(sum(v3.^2));

% --- Get XYZ coord of center of box ---
maxp = max(c);	% XYZmax vector
minp = min(c);	% XYZmin vector
origin = (maxp + minp)./2;
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function prescrip = origin_prescrip(origin)
% Return string with Siemens-style prescription for slice origin
%	e.q. L10 P5 H12.5 
%------------------------------------------------------------------------
c1 = 'L'; c2 = 'P'; c3 = 'H';
if (origin(1) < 0), c1 = 'R' ; end
if (origin(2) < 0), c2 = 'A' ; end
if (origin(3) < 0), c3 = 'F' ; end
prescrip = sprintf('%s%1.1f %s%1.1f %s%1.1f',c1,abs(origin(1)),c2,abs(origin(2)),c3,abs(origin(3)));
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
function plane = vector_prescrip(v3)
% Return slice orientation with Siemens-style prescription
%	e.q. T>C-20
% Inputs is slice normal vector.
%
% M.Elliott 4/16
%------------------------------------------------------------------------

% --- Get angles of slice normal to the X,Y,Z vectors ---
theta = [acos(v3' * [1 ; 0 ; 0]) ; acos(v3' * [0 ; 1 ; 0]) ;acos(v3' * [0 ; 0 ; 1])];
mag   = abs(cos(theta));
names = ['SAG' ; 'COR' ; 'TRA'];
inits = ['S' ; 'C' ; 'T'];
[~,order] = sort(mag);

% --- Find how many of the XYZ vectors are perpendicular to the image plane ---
tmp = find(mag < 0.01);
nperp = numel(tmp);
switch nperp
	case 0					% No perp vectors means DOUBLE oblique angle
		angle0 = rad2deg(theta(order(3)));

		plane = [inits(order(3)) '>' inits(order(2))];
		angle1 = rad2deg(theta(order(2))) - 90;
		if (angle0 > 90.) , angle1 = -angle1; end
		plane = sprintf('%s%1.1f',plane,angle1);

		plane = [plane '>' inits(order(1))];
		angle2 = rad2deg(theta(order(1))) - 90;
		if (angle0 > 90.) , angle2 = -angle2; end
		plane = sprintf('%s%1.1f',plane,angle2);
		
	case 1					% 1 perp vector means Single oblique angle
		plane = [inits(order(3)) '>' inits(order(2))];
		angle1 = rad2deg(theta(order(3)));
		angle2 = rad2deg(theta(order(2)));
		if (angle1 > 90.) , angle1 = angle1-180; end
		if (angle2 < 90.) , angle1 = -angle1; end
		plane = sprintf('%s%1.1f',plane,angle1);
		
	case 2					% 2 perp vectors means orthogonal slice (no oblique)
		plane = names(order(3),:);
end
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function res = nifti_grid(volinfo, include_rgrid, use_spmstyle)
% Returns coords of box around middle, top and bottom slice of a nifit volume
% Also returns (x,y,z) locations of all voxels in a Nifti Volume

if (nargin < 1) || isempty(include_rgrid), include_rgrid = 0; end
if (nargin < 2) || isempty(use_spmstyle),  use_spmstyle  = 1; end

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
if (include_rgrid)
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
    res.rvol = reshape(r,3,nx,ny,nz);
else
    res.rvol = [];
end

% return as a structure
res.cmid = c;
res.ctop = ctop;
res.cbot = cbot;
res.mat  = mat;
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function res = rda_3dcalc(info,xform)
% Return XYZ coords of bounding box of a SVS RDA file
% inputs:
%	info = the info structure from readsRDA2()
%
% outputs:
%   c    = [5,3] array of XYZ coords of the corners of the middle slice/plane of FOV (in magnet coords)
%           (really only 4 corners, 1st corner point is replicated as the 5th to facilitate calling plot3())
%   ctop = XYZ coords of the corners of the "top" slice of the FOV (or top plane if SVS)
%   cbot = XYZ coords of the corners of the "bootom" slice
%
% M.Elliott 06/2020
%------------------------------------------------------------------------
if (nargin < 2) || isempty(xform), xform = []; end

% get vectors for row, column and perp image directions
v1 = [info.RowVector0; info.RowVector1; info.RowVector2];
v2 = [info.ColumnVector0; info.ColumnVector1; info.ColumnVector2];
v3 = cross(v1,v2);

% normalize directions vectors in case they aren't already
v1 = v1/sqrt(sum(v1.^2));
v2 = v2/sqrt(sum(v2.^2));
v3 = v3/sqrt(sum(v3.^2));

pixelspacing = [info.PixelSpacingRow; info.PixelSpacingCol];
%columns      = info.NumberOfColumns;
%rows         = info.NumberOfRows;
columns      = checkstruct(info,'NumberOfColumns',1);   % VE11E RDA changes!!
rows         = checkstruct(info,'NumberOfRows',1);
xfov         = pixelspacing(1)*double(columns);
yfov         = pixelspacing(2)*double(rows);

% SVS are always rotated 90 degrees compared to images ??!!
tmp  = xfov;    % Actually it's just xfov/yfov are switched
xfov = yfov;
yfov = tmp;
 
% compute corners of bounding box, XYZ triples
c1 = [0 0 0]';
c2 = v1*xfov;
c3 = c2 + v2*yfov;
c4 = c1 + v2*yfov;
c = [c1' ; c2' ; c3' ; c4' ; c1'];	% add 1st point at end to close box

% translate box according to bottom corner info
corner = [info.PositionVector0; info.PositionVector1; info.PositionVector2];
c(:,1) = c(:,1) + corner(1);
c(:,2) = c(:,2) + corner(2);
c(:,3) = c(:,3) + corner(3);

% make boxes for top and bottom of slice
thick = info.SliceThickness;
ctop = c;
ctop(:,1) = ctop(:,1) + v3(1)*thick/2;
ctop(:,2) = ctop(:,2) + v3(2)*thick/2;
ctop(:,3) = ctop(:,3) + v3(3)*thick/2;
cbot = c;
cbot(:,1) = cbot(:,1) - v3(1)*thick/2;
cbot(:,2) = cbot(:,2) - v3(2)*thick/2;
cbot(:,3) = cbot(:,3) - v3(3)*thick/2;

if (~isempty(xform))
    c    = (xform*c')';    % put SVS bounding boxes in another space (probably NIFTI)
    ctop = (xform*ctop')';
    cbot = (xform*cbot')';
end

% return as a structure
res.cmid = c;
res.ctop = ctop;
res.cbot = cbot;
res.mat  = [];
res.rvol = [];
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------s
function res = dicom_3dcalc(info,mosaic_slice,mz,xform)
% Return XYZ coords of bounding box of a dicom image or SVS spectrum
% inputs:
%	info = the info structure from dicominfo()
%	mosaic_slice = -1, get box of first slice in mosaic
%				    0, get box of middle slice in mosaic
%				    1, get box of last slice in mosaic
%					(default = 0)
%   mz = number of slices in mosaic (if omitted, will be determined)
%
% outputs:
%   c    = [5,3] array of XYZ coords of the corners of the middle slice/plane of FOV (in magnet coords)
%           (really only 4 corners, 1st corner point is replicated as the 5th to facilitate calling plot3())
%   ctop = XYZ coords of the corners of the "top" slice of the FOV (or top plane if SVS)
%   cbot = XYZ coords of the corners of the "bootom" slice
%
% M.Elliott 2/07
%------------------------------------------------------------------------

% set whether we first, middle, or last slice in mosaic'd images
if (nargin < 2) || isempty(mosaic_slice), mosaic_slice = 0; end
if (nargin < 3) || isempty(mz),           mz           = -1; end
if (nargin < 4) || isempty(xform),       xform         = []; end

% --- is this an MRS file? ---
is_spec = isfield(info,'Private_7fe1_1010');

% get vectors for row, column and perp image directions
%%v1 = info.ImageOrientationPatient(1:3);
%%v2 = info.ImageOrientationPatient(4:6);
[stat,v,info] = dicom_get_header(info,'ImageOrientationPatient');	% This checks Siemens Private header 
v1 = v(1:3);
v2 = v(4:6);
v3 = cross(v1,v2);

% normalize directions vectors in case they aren't already
v1 = v1/sqrt(sum(v1.^2));
v2 = v2/sqrt(sum(v2.^2));
v3 = v3/sqrt(sum(v3.^2));

% fix info for Mosaic'd images
if (~isempty(strfind(info.ImageType,'MOSAIC')))
	nx = double(info.Width);					% Mosaic image size
	ny = double(info.Height);
	mx = double(info.AcquisitionMatrix(1));		% Actual image size
	my = double(info.AcquisitionMatrix(4));
	if (mx == 0)								% Handle strangeness??
		mx = double(info.AcquisitionMatrix(3));	
		my = double(info.AcquisitionMatrix(2));
	end
	if (mod(nx,mx) ~= 0)
		error('Dicom error in file:\n%s\n    Mosaiced X image dimension (%1d) is not an integer multiple of acq matrix size (%1d).\n',info.Filename,nx,mx);
	end
	if (mod(ny,my) ~= 0)
		error('Dicom error in file:\n%s\n    Mosaiced Y image dimension (%1d) is not an integer multiple of acq matrix size (%1d).\n',info.Filename,ny,my);
	end	
	% Find out how many actual images in Mosaic
	if (mz == -1)
		xtiles  = (nx/mx);						% # of tiles in mosaic in x-dir			
		ytiles  = (ny/my);
		% use thumbnail image to figure how many slices
		if (isfield(info,'IconImageSequence'))	
			[im,tx,ty] = dicom_thumb(info);		% get thumbnail image
			xtilsiz    = tx/xtiles;				% x-size of mosaic tiles
			ytilsiz    = ty/ytiles;
		% no thumbnail, have to use real image
		else
			im = dicomread(info);				% get real image
			xtilsiz    = mx;					% x-size of mosaic tiles
			ytilsiz    = my;			
		end
		bg      = min(im(:));					% background value (not always zero!)
		mz      = 0;
		for y=1:ytiles
			y0 = int16((y-0.5)*ytilsiz);
			for x = 1:xtiles
				x0 = int16((x-0.5)*xtilsiz);
%				fprintf(1,'%d %d %d\n',x0,y0,im(x0,y0));
				if (im(x0,y0) > bg), mz = mz+1; end
			end
		end
		fprintf(1,'Found %1d non-blank slices in mosaiced image.\n',mz);
	end
	dx = double(info.PixelSpacing(1));			% Voxel size
	dy = double(info.PixelSpacing(2));
	dz = double(info.SpacingBetweenSlices);		% (including slice gap)	
	info.ImagePositionPatient = info.ImagePositionPatient + v1 * (nx-mx)/2 * dx;
	info.ImagePositionPatient = info.ImagePositionPatient + v2 * (ny-my)/2 * dy;
	switch mosaic_slice
		case -1		% first slice
			info.ImagePositionPatient = info.ImagePositionPatient + v3 * 0 * dz;
		case  0		% middle slice
			info.ImagePositionPatient = info.ImagePositionPatient + v3 * (mz-1)/2 * dz;
		case  1		% last slice
			info.ImagePositionPatient = info.ImagePositionPatient + v3 * (mz-1) * dz;
	end
	xfov = dx * mx;
	yfov = dy * my;

% non-Mosaic: regular images or spectra
else
%%	xfov = double(info.PixelSpacing(1) * info.Width);
%%	yfov = double(info.PixelSpacing(2) * info.Height);
	[stat,pixelspacing,info] = dicom_get_header(info,'PixelSpacing');	% This also checks Siemens Private header (for Spectro files)
	[stat,columns,info]      = dicom_get_header(info,'Columns');		% same as "Width"
	[stat,rows,info]         = dicom_get_header(info,'Rows');			% same as "Hieght"
	xfov                     = pixelspacing(1)*double(columns);
	yfov                     = pixelspacing(2)*double(rows);

    % SVS are always rotated 90 degrees compared to images ??!!
    if(is_spec)
        tmp  = xfov;    % Actually it's just xfov/yfov are switched
        xfov = yfov;
        yfov = tmp;
    end
end

% compute corners of bounding box, XYZ triples
c1 = [0 0 0]';
c2 = v1*xfov;
c3 = c2 + v2*yfov;
c4 = c1 + v2*yfov;
c = [c1' ; c2' ; c3' ; c4' ; c1'];	% add 1st point at end to close box

% translate box according to bottom corner info
[stat,corner,info] = dicom_get_header(info,'ImagePositionPatient');	
%%c(:,1) = c(:,1) + info.ImagePositionPatient(1);
%%c(:,2) = c(:,2) + info.ImagePositionPatient(2);
%%c(:,3) = c(:,3) + info.ImagePositionPatient(3);
c(:,1) = c(:,1) + corner(1);
c(:,2) = c(:,2) + corner(2);
c(:,3) = c(:,3) + corner(3);

% make boxes for top and bottom of slice
[stat,thick,info] = dicom_get_header(info,'SpacingBetweenSlices');	
if (~stat || isempty(thick))
    [stat,thick,info] = dicom_get_header(info,'SliceThickness');	 
end
ctop = c;
ctop(:,1) = ctop(:,1) + v3(1)*thick/2;
ctop(:,2) = ctop(:,2) + v3(2)*thick/2;
ctop(:,3) = ctop(:,3) + v3(3)*thick/2;
cbot = c;
cbot(:,1) = cbot(:,1) - v3(1)*thick/2;
cbot(:,2) = cbot(:,2) - v3(2)*thick/2;
cbot(:,3) = cbot(:,3) - v3(3)*thick/2;

if (~isempty(xform))
    c    = (xform*c')';    % put SVS bounding boxes in another space (probably NIFTI)
    ctop = (xform*ctop')';
    cbot = (xform*cbot')';
end

% return as a structure
res.cmid = c;
res.ctop = ctop;
res.cbot = cbot;
res.mat  = [];
res.rvol = [];
end






