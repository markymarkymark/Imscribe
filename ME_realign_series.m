%_______________________________________________________________________
function P = ME_realign_series(P,flags)
% Realign a time series of 3D images to the first of the series.
% FORMAT P = realign_series(P,flags)
% P  - a vector of volumes (see spm_vol)
%-----------------------------------------------------------------------
% P(i).mat is modified to reflect the modified position of the image i.
% The scaling (and offset) parameters are also set to contain the
% optimum scaling required to match the images.
%
% Cut from spm_realign.m by MELLIOTT
%_______________________________________________________________________

if numel(P)<2, return; end;

skip = sqrt(sum(P(1).mat(1:3,1:3).^2)).^(-1)*flags.sep;
d    = P(1).dim(1:3);                                                                                                                        
lkp = flags.lkp;
rand('state',0); % want the results to be consistant.
if d(3) < 3,
	lkp = [1 2 6];
	[x1,x2,x3] = ndgrid(1:skip(1):d(1)-.5, 1:skip(2):d(2)-.5, 1:skip(3):d(3));
	x1   = x1 + rand(size(x1))*0.5;
	x2   = x2 + rand(size(x2))*0.5;
else
	[x1,x2,x3]=ndgrid(1:skip(1):d(1)-.5, 1:skip(2):d(2)-.5, 1:skip(3):d(3)-.5);
	x1   = x1 + rand(size(x1))*0.5;
	x2   = x2 + rand(size(x2))*0.5;
	x3   = x3 + rand(size(x3))*0.5; 
end;

x1   = x1(:);
x2   = x2(:);
x3   = x3(:);

% Possibly mask an area of the sample volume.
%-----------------------------------------------------------------------
if ~isempty(flags.PW),
	[y1,y2,y3]=coords([0 0 0  0 0 0],P(1).mat,flags.PW.mat,x1,x2,x3);
	wt  = spm_sample_vol(flags.PW,y1,y2,y3,1);
	msk = find(wt>0.01);
	x1  = x1(msk);
	x2  = x2(msk);
	x3  = x3(msk);
	wt  = wt(msk);
else
	wt = [];
end;

% Compute rate of change of chi2 w.r.t changes in parameters (matrix A)
%-----------------------------------------------------------------------
V   = smooth_vol(P(1),flags.interp,flags.wrap,flags.fwhm);
deg = [flags.interp*[1 1 1]' flags.wrap(:)];

[G,dG1,dG2,dG3] = spm_bsplins(V,x1,x2,x3,deg);
clear V
A0 = make_A(P(1).mat,x1,x2,x3,dG1,dG2,dG3,wt,lkp);

b  = G;
if ~isempty(wt), b = b.*wt; end;

%-----------------------------------------------------------------------
if numel(P) > 2,
	% Remove voxels that contribute very little to the final estimate.
	% Simulated annealing or something similar could be used to
	% eliminate a better choice of voxels - but this way will do for
	% now. It basically involves removing the voxels that contribute
	% least to the determinant of the inverse covariance matrix.

	spm_chi2_plot('Init','Eliminating Unimportant Voxels',...
		      'Relative quality','Iteration');
	Alpha = spm_atranspa([A0 b]);
	det0  = det(Alpha);
	det1  = det0;
	spm_chi2_plot('Set',det1/det0);
	while det1/det0 > flags.quality,
		dets  = zeros(size(A0,1),1);
		for i=1:size(A0,1),
			dets(i) = det(Alpha - spm_atranspa([A0(i,:) b(i)]));
		end;
		[junk,msk] = sort(det1-dets);
		msk        = msk(1:round(length(dets)/10));
		 A0(msk,:) = [];   b(msk,:) = [];   G(msk,:) = [];
		 x1(msk,:) = [];  x2(msk,:) = [];  x3(msk,:) = [];
		dG1(msk,:) = []; dG2(msk,:) = []; dG3(msk,:) = [];
		if ~isempty(wt),  wt(msk,:) = []; end;
		Alpha = spm_atranspa([A0 b]);
		det1  = det(Alpha);
		spm_chi2_plot('Set',single(det1/det0));
	end;
	spm_chi2_plot('Clear');
end;
%-----------------------------------------------------------------------


if flags.rtm,
	count = ones(size(b));
	ave   = G;
	grad1 = dG1;
	grad2 = dG2;
	grad3 = dG3;
end;

spm_progress_bar('Init',length(P)-1,'Registering Images');
% Loop over images
%-----------------------------------------------------------------------
for i=2:length(P),
	V  = smooth_vol(P(i),flags.interp,flags.wrap,flags.fwhm);
	d  = [size(V) 1 1];
	d  = d(1:3);
	ss = Inf;
	countdown = -1;
	for iter=1:64,
		[y1,y2,y3] = coords([0 0 0  0 0 0],P(1).mat,P(i).mat,x1,x2,x3);
		msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
		if length(msk)<32, error_message(P(i)); end;

		F          = spm_bsplins(V, y1(msk),y2(msk),y3(msk),deg);
		if ~isempty(wt), F = F.*wt(msk); end;

		A          = A0(msk,:);
		b1         = b(msk);
		sc         = sum(b1)/sum(F);
		b1         = b1-F*sc;
%		soln       = spm_atranspa(A)\(A'*b1); % MELLIOTT
		soln       = (A'*A)\(A'*b1);

		p          = [0 0 0  0 0 0  1 1 1  0 0 0];
		p(lkp)     = p(lkp) + soln';
		P(i).mat   = inv(spm_matrix(p))*P(i).mat;

		pss        = ss;
		ss         = sum(b1.^2)/length(b1);
		if (pss-ss)/pss < 1e-8 && countdown == -1, % Stopped converging.
			countdown = 2;
		end;
		if countdown ~= -1,
			if countdown==0, break; end;
			countdown = countdown -1;
		end;
	end;
	if flags.rtm,
		% Generate mean and derivatives of mean
		tiny = 5e-2; % From spm_vol_utils.c
		msk        = find((y1>=(1-tiny) & y1<=(d(1)+tiny) &...
		                   y2>=(1-tiny) & y2<=(d(2)+tiny) &...
		                   y3>=(1-tiny) & y3<=(d(3)+tiny)));
		count(msk) = count(msk) + 1;
		[G,dG1,dG2,dG3] = spm_bsplins(V,y1(msk),y2(msk),y3(msk),deg);
		ave(msk)   = ave(msk)   +   G*sc;
		grad1(msk) = grad1(msk) + dG1*sc;
		grad2(msk) = grad2(msk) + dG2*sc;
		grad3(msk) = grad3(msk) + dG3*sc;
	end;
	spm_progress_bar('Set',i-1);
end;
spm_progress_bar('Clear');

if ~flags.rtm, return; end;
%_______________________________________________________________________
M=P(1).mat;
A0 = make_A(M,x1,x2,x3,grad1./count,grad2./count,grad3./count,wt,lkp);
if ~isempty(wt), b = (ave./count).*wt;
else b = (ave./count); end

clear ave grad1 grad2 grad3

% Loop over images
%-----------------------------------------------------------------------
spm_progress_bar('Init',length(P),'Registering Images to Mean');
for i=1:length(P),
	V  = smooth_vol(P(i),flags.interp,flags.wrap,flags.fwhm);
	d  = [size(V) 1 1 1];
	ss = Inf;
	countdown = -1;
	for iter=1:64,
		[y1,y2,y3] = coords([0 0 0  0 0 0],M,P(i).mat,x1,x2,x3);
		msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
		if length(msk)<32, error_message(P(i)); end;

		F          = spm_bsplins(V, y1(msk),y2(msk),y3(msk),deg);
		if ~isempty(wt), F = F.*wt(msk); end;

		A          = A0(msk,:);
		b1         = b(msk);
		sc         = sum(b1)/sum(F);
		b1         = b1-F*sc;
		soln       = spm_atranspa(A)\(A'*b1);

		p          = [0 0 0  0 0 0  1 1 1  0 0 0];
		p(lkp)     = p(lkp) + soln';
		P(i).mat   = inv(spm_matrix(p))*P(i).mat;

		pss        = ss;
		ss         = sum(b1.^2)/length(b1);
		if (pss-ss)/pss < 1e-8 && countdown == -1 % Stopped converging.
			% Do three final iterations to finish off with
			countdown = 2;
		end;
		if countdown ~= -1
			if countdown==0, break; end;
			countdown = countdown -1;
		end;
	end;
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');


% Since we are supposed to be aligning everything to the first
% image, then we had better do so
%-----------------------------------------------------------------------
M = M/P(1).mat;
for i=1:length(P)
	P(i).mat   = M*P(i).mat;
end

return;
%_______________________________________________________________________

%_______________________________________________________________________
function [y1,y2,y3]=coords(p,M1,M2,x1,x2,x3)
% Rigid body transformation of a set of coordinates.
M  = (inv(M2)*inv(spm_matrix(p))*M1);
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function V = smooth_vol(P,hld,wrp,fwhm)
% Convolve the volume in memory.
s  = sqrt(sum(P.mat(1:3,1:3).^2)).^(-1)*(fwhm/sqrt(8*log(2)));
x  = round(6*s(1)); x = -x:x;
y  = round(6*s(2)); y = -y:y;
z  = round(6*s(3)); z = -z:z;
x  = exp(-(x).^2/(2*(s(1)).^2));
y  = exp(-(y).^2/(2*(s(2)).^2));
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
d  = [hld*[1 1 1]' wrp(:)];
V  = spm_bsplinc(P,d);
spm_conv_vol(V,V,x,y,z,-[i j k]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function A = make_A(M,x1,x2,x3,dG1,dG2,dG3,wt,lkp)
% Matrix of rate of change of weighted difference w.r.t. parameter changes
p0 = [0 0 0  0 0 0  1 1 1  0 0 0];
A  = zeros(numel(x1),length(lkp));
for i=1:length(lkp)
	pt         = p0;
	pt(lkp(i)) = pt(i)+1e-6;
	[y1,y2,y3] = coords(pt,M,M,x1,x2,x3);
	tmp        = sum([y1-x1 y2-x2 y3-x3].*[dG1 dG2 dG3],2)/(-1e-6);
	if ~isempty(wt), A(:,i) = tmp.*wt;
	else A(:,i) = tmp; end
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function error_message(P)
str = {	'There is not enough overlap in the images',...
	'to obtain a solution.',...
	' ',...
	'Offending image:',...
	 P.fname,...
	' ',...
	'Please check that your header information is OK.'};
spm('alert*',str,mfilename,sqrt(-1));
error('insufficient image overlap')
%_______________________________________________________________________

