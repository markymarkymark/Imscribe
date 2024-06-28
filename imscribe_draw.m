function imscribe_draw(VOL1,VOL2,VOL3,c1,c2,c3,cx,opt,winindex)
% -------------------------------------------------------------------------
% 1) In winindex(1): Draw VOL1 w/ VOL2 contour overlay
% 2) In winindex(2): Draw VOL3 w/ cx FOV ovverlay
% -------------------------------------------------------------------------

% Init color map
cmap           = colormap(gray(256));
nc             = size(cmap,1);
contour_levels = round(opt.contourlevel/100*10)+1;
opt.txtid(winindex(3)).Value = cx.prescrip;

% --- WININDEX(1) ---------------------------------
% --- Show the VOL1 w/ countour overlay of VOL2 ---
axes(opt.winid(winindex(1))); 
I1 = extractslice(VOL1,c1.rvol,opt.slicenum,opt.imbright,c1.prescrip(1),opt.displane);
imshow(flipim(I1,c1.prescrip(1),opt.displane),[0 1]);
if (~isempty(VOL2))
    I2 = extractslice(VOL2,c2.rvol,opt.slicenum,opt.imbright,c2.prescrip(1),opt.displane);
    hold on;
    contour(flipim(I2,c2.prescrip(1),opt.displane),contour_levels,'r');
    hold off;
    drawnow;
end

% --- WININDEX(2) ----------------------
% --- Show VOL3 w/ FOV of CX overlay --- 
axes(opt.winid(winindex(2)));
[I3,~,xv3,yv3,zv3] = extractslice(VOL3,c3.rvol,opt.slicenum,opt.imbright,c3.prescrip(1),opt.displane);
if (min(size(xv3)) < 2), fprintf(2,'WARNING: Cannot replane single slice Structural. Skipping display.\n'); return; end
me_warp(xv3,yv3,zv3,I3); 
hold on
plot_box(cx.ctop,cx.cbot,'g',2);
plot_box(cx.cmid,cx.cmid,'g',1);  % middle slice of FOV w/ thin line
hfill = fill3(cx.cmid(1:4,1),cx.cmid(1:4,2),cx.cmid(1:4,3),'g');
hfill.FaceAlpha=0.2;

% --- Specify viewing angle ---
setview(opt.displane);
hold off;
axis off;
%axis([xr yr zr]*1.1);
drawnow;
h1 = rotate3d(opt.winid(winindex(2)));    % enable rotating of view with mouse
set(h1,'Enable','on');
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function plot_box(ctop,cbot,color,width)
if (nargin < 3), color = 'k'; end
if (nargin < 4), widhth = 1; end

plot3(cbot(:,1),cbot(:,2),cbot(:,3), color,'LineWidth',width);
plot3(ctop(:,1),ctop(:,2),ctop(:,3), color,'LineWidth',width);
for j=1:4, plot3([cbot(j,1);ctop(j,1)], [cbot(j,2);ctop(j,2)], [cbot(j,3);ctop(j,3)],color,'LineWidth',width); end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function imout = flipim(imin,orient_code,plane)
% Flip a 2D image so it displays correctly w/ imshow
switch orient_code      % orientation of 3D volume the 2d image came from
    case 'T'
        switch plane
            case 1      % Axial slice wanted
                imout = flipdim(imin',1);
            case 2      % Coronal slice wanted
                imout = flipdim(imin',1);
            case 3      % Sagittal slice wanted
                imout = flipdim(flipdim(imin',1),2);
        end
        
    case 'C'
        switch plane
            case 1      % Axial slice wanted
                imout = imin';
            case 2      % Coronal slice wanted
                imout = flipdim(imin',1);
            case 3      % Sagittal slice wanted
                imout = flipdim(imin,1);
        end
        
    case 'S'
        switch plane
            case 1      % Axial slice wanted
                imout = flipdim(imin,2);
            case 2      % Coronal slice wanted
                imout = flipdim(flipdim(imin,1),2);
            case 3      % Sagittal slice wanted
                imout = flipdim(imin',1);
        end
end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [im,s0,xv,yv,zv] = extractslice(vol,r,slice,bright,orient_code,plane)
% Pull slice from 3D volume. Get desired plane given orientation of 3D input vol
%   slice is expressed in % of number of slices in volume
%
%   if slice = [i,j,k], then use index correct for desired orientation of display
[n1,n2,n3] = size(vol);

switch orient_code      % orientation of 3D volume
    case 'T'
        switch plane
            case 1      % Axial slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n3-1)) + 1; % use slice = % of FOV
                else                    s0 = slice(3);                    % use k from [i,j,k]
                end
                s0 = min([max([s0 1]) n3]);
                im = vol(:,:,s0);
                if (nargout > 2)
                    xv = squeeze(r(1,:,:,s0));
                    yv = squeeze(r(2,:,:,s0));
                    zv = squeeze(r(3,:,:,s0));
                end
            case 2      % Coronal slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n2-1)) + 1; % use slice = % of FOV
                else                    s0 = slice(2);                    % use j from [i,j,k]
                end
                s0 = min([max([s0 1]) n2]);
                im = squeeze(vol(:,s0,:));
                if (nargout > 2)
                    xv = squeeze(r(1,:,s0,:));
                    yv = squeeze(r(2,:,s0,:));
                    zv = squeeze(r(3,:,s0,:));
                end
            case 3      % Sagittal slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n1-1)) + 1; % use slice = % of FOV
                else                    s0 = slice(1);                    % use i from [i,j,k]
                end
                s0 = min([max([s0 1]) n1]);
                im = squeeze(vol(s0,:,:));
                if (nargout > 2)
                    xv = squeeze(r(1,s0,:,:));
                    yv = squeeze(r(2,s0,:,:));
                    zv = squeeze(r(3,s0,:,:));
                end
        end
        
    case 'C'
        switch plane
            case 1      % Axial slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n2-1)) + 1;
                else                    s0 = slice(2);
                end
                im = squeeze(vol(:,s0,:));
                if (nargout > 2)
                    xv = squeeze(r(1,:,s0,:));
                    yv = squeeze(r(2,:,s0,:));
                    zv = squeeze(r(3,:,s0,:));
                end
            case 2      % Coronal slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n3-1)) + 1;
                else                    s0 = slice(3);
                end
                im = vol(:,:,s0);
                if (nargout > 2)
                    xv = squeeze(r(1,:,:,s0));
                    yv = squeeze(r(2,:,:,s0));
                    zv = squeeze(r(3,:,:,s0));
                end
            case 3      % Sagittal slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n1-1)) + 1;
                else                    s0 = slice(1);
                end
                im = squeeze(vol(s0,:,:));
                if (nargout > 2)
                    xv = squeeze(r(1,s0,:,:));
                    yv = squeeze(r(2,s0,:,:));
                    zv = squeeze(r(3,s0,:,:));
                end
        end
        
    case 'S'
        switch plane
            case 1      % Axial slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n2-1)) + 1; % use slice = % of FOV
                else                    s0 = slice(2);                    % use k from [i,j,k]
                end
                im = squeeze(vol(:,s0,:));
                if (nargout > 2)
                    xv = squeeze(r(1,:,s0,:));
                    yv = squeeze(r(2,:,s0,:));
                    zv = squeeze(r(3,:,s0,:));
                end
            case 2      % Coronal slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n1-1)) + 1; % use slice = % of FOV
                else                    s0 = slice(1);                    % use i from [i,j,k]
                end
                im = squeeze(vol(s0,:,:));
                if (nargout > 2)
                    xv = squeeze(r(1,s0,:,:));
                    yv = squeeze(r(2,s0,:,:));
                    zv = squeeze(r(3,s0,:,:));
                end
            case 3      % Sagittal slice wanted
                if (numel(slice) == 1), s0 = round(slice/100*(n3-1)) + 1; % use slice = % of FOV
                else                    s0 = slice(3);                    % use k from [i,j,k]
                end
                im = vol(:,:,s0);
                if (nargout > 2)
                    xv = squeeze(r(1,:,:,s0));
                    yv = squeeze(r(2,:,:,s0));
                    zv = squeeze(r(3,:,:,s0));
                end
        end
end

% brighten/darken image using "bright" factor
im = im/max(im(:))*bright;
if (bright > 1)
    im = clip(im,0,1);
end
end

% -------------------------------------------------------------------------
% clip input array to min,max limits
%  (uses matrix Logicals to replace Loops)
function  y = clip( x, lo, hi )
y = (x .* [x<=hi])  +  (hi .* [x>hi]);
y = (y .* [x>=lo])  +  (lo .* [x<lo]);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function angles = setview(plane)
switch plane
    case 1      % Axial slice
        angles=[180,90];
    case 2      % Coronal slice
        angles=[180,0];
    case 3      % Sagittal slice
        angles=[89.9,0.1];
end
view(angles);
end
