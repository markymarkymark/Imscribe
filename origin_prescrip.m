function prescrip = origin_prescrip(origin)
% Return string with Siemens-style prescription for slice origin
%	e.q. L10 P5 H12.5 
%------------------------------------------------------------------------
c1 = 'L'; c2 = 'P'; c3 = 'H';
if (origin(1) < 0), c1 = 'R' ; end
if (origin(2) < 0), c2 = 'A' ; end
if (origin(3) < 0), c3 = 'F' ; end
prescrip = sprintf('%s%1.1f %s%1.1f %s%1.1f',c1,abs(origin(1)),c2,abs(origin(2)),c3,abs(origin(3)));
return