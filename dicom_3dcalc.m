function [c,ctop,cbot] = dicom_3dcalc(info,mosaic_slice,mz)
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
if nargin < 2 , mosaic_slice = 0; end
if nargin < 3 , mz = -1; end

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

end



