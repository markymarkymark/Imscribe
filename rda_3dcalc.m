function [c,ctop,cbot] = rda_3dcalc(info)
% Return XYZ coords of bounding box of a Siemens SVS RDA file
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

end



