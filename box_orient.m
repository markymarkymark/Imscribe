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

return
