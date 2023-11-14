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

return
