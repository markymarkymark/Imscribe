function [prescrip,origin,rot] = box_prescrip(c)
% Return string with Siemens-style prescription
%	e.q. L10 P5 H12.5 T>C-20
% Input is bounding box coords, as returned by dicom_3dcalc(), or nifti_grid()
%   NOTE: the bounding box should be in Siemens coord space (i.e. nifti_grid() result needs mat_N2S*)
%
% M.Elliott 2/07
%------------------------------------------------------------------------

% --- get inplane & normal vectors, and origin of bounding box ---
[v1,v2,v3,origin] = box_orient(c);

% --- Get Siemens desciption of slice plane orientation ---
plane = vector_prescrip(v3);

% --- Add center voxel location ---
center = origin_prescrip(origin);
% c1 = 'L'; c2 = 'P'; c3 = 'H';
% if (origin(1) < 0), c1 = 'R' ; end
% if (origin(2) < 0), c2 = 'A' ; end
% if (origin(3) < 0), c3 = 'F' ; end
% prescrip = sprintf('%s %s%1.1f %s%1.1f %s%1.1f',plane,c1,abs(origin(1)),c2,abs(origin(2)),c3,abs(origin(3)));
prescrip = [plane ' ' center];

% --- Determine the inplane rotation ---
rot = inplane_rot(prescrip,v1,v2,v3);
return

