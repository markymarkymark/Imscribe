function [plane] = vector_prescrip(v3)
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

return

