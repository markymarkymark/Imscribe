function cout = transform_coords(cin,mat)
% convert matrix of XYZ coords to new space mapped to by Affine mat

ctmp = [cin ones(5,1)];     % make each XYZ into XYZ1
cout = mat * ctmp';			% transform to TARGET space
cout = cout(1:3,:)';		% clip XYZ1 into XYZ
return;

