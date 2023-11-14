% imscribe_engine()
function [status,p2x] = imscribe_engine(paths,files1,files2,files3,options,winid,texthandle)

% --- Enable re-use of previous xform if desired ---
global last_Mat last_rvol3
global vol1 vol2 vol3 is_spec is_rda c2 ctop2 cbot2   % for re-use when redraw option is on

% --- Set return variables ---
status = 0;	% this means error
p2x    = '';

% --- check args ---
if (nargin < 6), winid = []; end
if (nargin < 7), texthandle = []; end

% --- GUI - clear previous contents of Siemens prescriptions ---
if (~isempty(winid))
    set(texthandle(1),'String','');
    set(texthandle(2),'String','');
    
    % -- reset plot windows ---
    cla(winid(1)); title(winid(1),''); axis(winid(1),'off');
    cla(winid(2)); title(winid(2),''); axis(winid(2),'off');
    cla(winid(3)); title(winid(3),''); axis(winid(3),'off');
    cla(winid(4)); title(winid(4),''); axis(winid(4),'off');
    
    h1 = rotate3d(winid(1));
    h3 = rotate3d(winid(3));
    set(h1,'Enable','off');
    set(h3,'Enable','off');
    drawnow
    
    cmap           = colormap(gray(256));
    nc             = size(cmap,1);
    contour_levels = round(options.contourlevel/100*10)+1;
end

% --- Get paths for possible output images ---
scratchpath = [paths{1} filesep()]; % scratch folder for temp Nifti files
svsmaskpath = [paths{2} filesep()]; % location for Nifti with mask of SVS ROI
dicompath   = [paths{3} filesep()]; % location for Dicoms with ROI mask

% --- Set some local variables ---
mat_N2S        = [-1 0 0; 0 -1 0 ; 0 0 1]'; % xform from Nifti to Siemens MRI coord space ---
mat_S2N        = inv(mat_N2S);

% --- Check for illegal combinations of options ---
if ((options.reuse_xform) && (isempty(last_Mat)))
    fprintf(2,'ERROR: No previous xformation matrix exists to re-use.\n');
    return;
end
if (options.redraw)
    if (options.svsmask && isempty(vol1))
        fprintf(2,'ERROR: Missing previous Template data to re-draw!\n');
        return;
    end
    if (~options.svsmask && (isempty(vol1) || isempty(vol3)))
        fprintf(2,'ERROR: Missing previous Template and Target data to re-draw!\n');
        return;
    end
end
if (options.dicom_mask && options.is_nifti(3))
    fprintf(2,'ERROR: Can only generate masked Dicoms if TARGET localizer is from Dicom files.\n');
    return;
end
if (options.dicom_mask && ~options.is_ROI)
    fprintf(2,'ERROR: Will only generate masked Dicoms if "select ROI" if on.\n');
    return;
end

% --- Start counting program duration ---
tic;

disp(' ');
disp('---------------------------------------------------------------------------------');
disp('Starting ImScribe...');
fprintf(1,'   TEMPLATE Localizer: %s\n',files1{1});
fprintf(1,'   TEMPLATE FOV/ROI  : %s\n',files2{1});
if (~options.svsmask && ~options.template_only)
    fprintf(1,'   TARGET Localizer:   %s\n',files3{1});
end

% --- Convert Dicoms to Nifti with SPM, or just read Niftis ---
if (~options.redraw)
    disp('---------------------------------------------------------------------------------');
    disp('Reading TEMPLATE and TARGET files...');
    
    %--- TEMPLATE LOCALIZER ---
    niftifile1 = [scratchpath 'TemplateLoc.nii'];
    if (~options.is_nifti(1))
        fprintf(1,'   Converting %1d TEMPLATE Localizer Dicoms to Nifti\n',numel(files1));
        vol1       = dicom_spm2nifti(files1,niftifile1);
        if (iscell(vol1)), vol1 = vol1{1}; end
    else
        fprintf(1,'   Reading TEMPLATE Localizer Nifti header\n');
        vol1 = spm_vol(files1{1});
        [stat,message,messageid] = copyfile(files1{1},niftifile1,'f');
        fprintf(1,'   Saved a copy of the TEMPLATE Localizer to %s\n',niftifile1);
    end
    
    % --- Check if TEMPLATE FOV Dicom file(s) are RDA file, or Spectro, or Mosaic'd ---
    if (options.is_rda)
        is_rda    = 1;
        is_spec   = 1;
        is_mosaic = 0;
        [~,rdahdr] = readRDA2(files2{1});
    elseif (~options.is_nifti(2))
        is_rda    = 0;
        dcmhdr    = dicominfo(files2{1},'UseDictionaryVR',true);
        is_spec   = isfield(dcmhdr,'Private_7fe1_1010') || isfield(dcmhdr,'SpectroscopyData');
        is_mosaic = ~isempty(strfind(dcmhdr.ImageType,'MOSAIC'));
        if ((is_spec || is_mosaic) && numel(files2) > 1)
            fprintf(2,'   WARNING: Multiple Spectro or Mosaic Dicom files selected for TEMPLATE FOV. Using first file only!\n');
            tmp = cell(1); tmp{1} = files2{1}; files2 = tmp;
        end
    else
        is_spec   = 0;
        is_mosaic = 0;
    end
    
    % --- SVSmask option requires dicom or RDA from SVS data ---
    if (options.svsmask && ~is_spec)
        fprintf(2,'ERROR: SVS mask option requires TEMPLATE FOV from an SVS dicom or RDA file.\n');
        return;
    end
    
    %--- TEMPLATE FOV or ROI ---
    if (options.is_ROI), niftifile2 = [scratchpath 'TemplateROI.nii'];
    else                 niftifile2 = [scratchpath 'TemplateFOV.nii']; end
    if (~options.is_nifti(2))
        
        % --- Convert Image Dicom(s) to Nifti ---
        if (~is_spec)
            fprintf(1,'   Converting %1d TEMPLATE FOV Dicoms to Nifti\n',numel(files2));
            niftifile2 = [scratchpath 'TemplateFOV.nii'];
            vol2       = dicom_spm2nifti(files2,niftifile2);
            if (iscell(vol2)), vol2 = vol2{1}; end
            
        % --- Find Dicom or RDA SVS volume location ---
        else
            if (options.is_ROI)
                fprintf(2,'   WARNING: ROI mode is incompatible with Spectro TEMPLATE. Ignoring ROI setting.\n');
                options.is_ROI = 0;
            end
            if (~is_rda)
                fprintf(1,'   Calculating Dicom SVS location\n');
                [c2,ctop2,cbot2] = dicom_3dcalc(dcmhdr);
                c2    = (mat_S2N*c2')';    % put bounding boxes in Nifti space
                ctop2 = (mat_S2N*ctop2')';
                cbot2 = (mat_S2N*cbot2')';
            else % RDA file!!!
                fprintf(1,'   Calculating RDA file SVS location\n');
                [c2,ctop2,cbot2] = rda_3dcalc(rdahdr);
                c2    = (mat_S2N*c2')';    % put bounding boxes in Nifti space
                ctop2 = (mat_S2N*ctop2')';
                cbot2 = (mat_S2N*cbot2')';
            end                
        end
    else
        fprintf(1,'   Reading TEMPLATE FOV/ROI Nifti header\n');
        vol2 = spm_vol(files2{1});
        [stat,message,messageid] = copyfile(files2{1},niftifile2,'f');
        fprintf(1,'   Saved a copy of the TEMPLATE FOV/ROI to %s\n',niftifile2);
    end
    if (numel(vol2) > 1), vol2 = vol2(1); end   % in case 4D Nifti was chosen
    
    %--- TARGET LOCALIZER ---
    niftifile3 = [scratchpath 'TargetLoc.nii'];
    if (~options.svsmask && ~options.template_only)
        if (~options.is_nifti(3))
            fprintf(1,'   Converting %1d TARGET Localizer Dicoms to Nifti\n',numel(files3));
            vol3 = dicom_spm2nifti(files3,niftifile3);
            if (iscell(vol3)), vol3 = vol3{1}; end
        else
            fprintf(1,'   Reading TARGET Localizer Nifti header\n');
            vol3 = spm_vol(files3{1});
            [stat,message,messageid] = copyfile(files3{1},niftifile3,'f');
            fprintf(1,'   Saved a copy of the TARGET Localizer to %s\n',niftifile3);
        end
    end
    
    % --- Skull strip Localizers ---
    if (options.skullstrip_temploc)
       	fprintf(1,'   Skull stripping Template Localizer...\n');
        [~,ssroot] = fileparts(vol1.fname);
        ssroot     = [ssroot '_SS'];
        outfile    = ME_spm8_skullstrip(vol1.fname,ssroot);
       	fprintf(1,'   ... saved result to %s\n\n',outfile);
        vol1       = spm_vol(outfile);
    end
    if (options.skullstrip_targloc)
       	fprintf(1,'   Skull stripping Target Localizer...\n');
        [~,ssroot] = fileparts(vol3.fname);
        ssroot     = [ssroot '_SS'];
        outfile    = ME_spm8_skullstrip(vol3.fname,ssroot);
       	fprintf(1,'   ... saved result to %s\n\n',outfile);
        vol3       = spm_vol(outfile);
    end
    
end % end re-draw skip

% --- Check that Localizers are at least 2 slices ---
% if (vol1.dim(3) < 2 || vol3.dim(3) < 2)
%     single_slice = 1;
%     %     fprintf(2,'ERROR: Don''t know how to coregister single slice Template or Target Localizers!\n');
%     %     return;
% else
%     single_slice = 0;
% end

% --- Compute Siemens-style "prescription" of TEMPLATES and TARGET images ---
disp('---------------------------------------------------------------------------------');
disp('Siemens-style Prescriptions...');

% --- Get grids & bounding boxes for the 3 volumes ---
[c1,ctop1,cbot1,r1vol] = nifti_grid(vol1);
if (~is_spec), [c2,ctop2,cbot2,r2vol] = nifti_grid(vol2); end
if (~options.svsmask && ~options.template_only), [c3,ctop3,cbot3,r3vol] = nifti_grid(vol3); end

% --- Get Siemens prescription from bounding box ---
[p1,o1,rot1] = box_prescrip((mat_N2S*c1')'); % convert to Siemens coord space
[p2,o2,rot2] = box_prescrip((mat_N2S*c2')');
if (~options.svsmask && ~options.template_only), [p3,o3,rot3] = box_prescrip((mat_N2S*c3')'); end
[lx,ly,lz]   = box_size((mat_N2S*c2')', (mat_N2S*ctop2')', (mat_N2S*cbot2')');
fprintf(1,'   TEMPLATE loc: %s \n',p1);
fprintf(1,'   TEMPLATE fov: %s (extent = %1.1f x %1.1f x %1.1f)\n',p2,lx,ly,lz);
if (~options.svsmask && ~options.template_only), fprintf(1,'   TARGET   loc: %s \n',p3); end

% --- Read in the voxel data ---
I1vol = spm_read_vols(vol1);
if (options.is_ROI), I2vol = spm_read_vols(vol2);	end % don't actually need voxel data from images or spectro
if (~options.svsmask && ~options.template_only), I3vol = spm_read_vols(vol3); end

% --- Handle Template FOV being an ROI (from Nifti file only) ---
if (options.is_ROI)
    disp('---------------------------------------------------------------------------------');
    disp('Finding TEMPLATE ROI extent...');
    
    % --- Find voxels in ROI ---
    if (options.ROIthresh > 0)
        fprintf(1,'   Finding voxels in TEMPLATE ROI that match value of %1d.\n',options.ROIthresh);
        proi = find(round(I2vol) == options.ROIthresh);
    else
        fprintf(1,'   Finding voxels in TEMPLATE ROI that are >= %1d.\n',abs(options.ROIthresh));
        proi = find((I2vol) >= abs(options.ROIthresh)); % negative value means threshold at >= |value|
    end
    roicount = numel(proi);
    if (roicount == 0)
        fprintf(2,'ERROR: No VALID voxels in Template ROI\n');
        return;
    end
    fprintf(1,'   Found %1d voxels in ROI.\n',roicount);
    
    % --- get slice through middle of ROI (for display) ---
    %[iroi,jroi,kroi]=ind2sub(size(I2vol),proi);              % decompose 1D index into (i,j,k)
    %mroi = round([mean(iroi); mean(jroi); mean(kroi)]);      % (i,j,k) center of ROI mass
    
    % --- get center of ROI and it's extent in XYZ ---
    rroi = [mean(r2vol(1,proi(:))); mean(r2vol(2,proi(:))); mean(r2vol(3,proi(:)))]; % (x,y,z) center of ROI mass
    xroi = [min(r2vol(1,proi(:))) max(r2vol(1,proi(:)))];                            % extent of ROI
    yroi = [min(r2vol(2,proi(:))) max(r2vol(2,proi(:)))];
    zroi = [min(r2vol(3,proi(:))) max(r2vol(3,proi(:)))];
    
    % --- get i,j,k location of ROI center in vol1 ---
    ijk1roi = round(inv(vol1.mat) * [rroi' 1]');  % i,j,k index into vol1 for location rroi
    
    % --- make box that encloses xyz extent of ROI ---
    b4 = [xroi(1) yroi(1) (zroi(2)+zroi(1))/2]';
    b3 = [xroi(2) yroi(1) (zroi(2)+zroi(1))/2]';
    b2 = [xroi(2) yroi(2) (zroi(2)+zroi(1))/2]';
    b1 = [xroi(1) yroi(2) (zroi(2)+zroi(1))/2]';
    c2roi    = [b1' ; b2' ; b3' ; b4' ; b1'];	% add 1st point at end to close box
    cbot2roi = c2roi; cbot2roi(:,3) = zroi(1);
    ctop2roi = c2roi; ctop2roi(:,3) = zroi(2);
    
    % --- find Siemens prescription for ROI ---
    [p2roi,o2roi,rot2roi] = box_prescrip((mat_N2S*c2roi')');
    [lx,ly,lz]            = box_size((mat_N2S*c2roi')', (mat_N2S*ctop2roi')', (mat_N2S*cbot2roi')');
    
    o2roix = mat_N2S*rroi;
    fprintf(1,'   TEMPLATE roi: %s\n',p2roi);
    fprintf(1,'   center of mass:   L%1.1f P%1.1f H%1.1f\n',o2roix(1),o2roix(2),o2roix(3));
    fprintf(1,'   extent:           %1.1f x %1.1f x %1.1f\n',lx,ly,lz);
    
    % --- find extent & center of SVS FOV ---
elseif (is_spec)
    svsx0 = (max([cbot2(:,1) ;ctop2(:,1)]) + min([cbot2(:,1) ;ctop2(:,1)]))/2;
    svsy0 = (max([cbot2(:,2) ;ctop2(:,2)]) + min([cbot2(:,2) ;ctop2(:,2)]))/2;
    svsz0 = (max([cbot2(:,3) ;ctop2(:,3)]) + min([cbot2(:,3) ;ctop2(:,3)]))/2;
    rsvs = [svsx0 ; svsy0 ; svsz0]; % (x,y,z) center of SVS FOV
    ijk1svs = round(inv(vol1.mat) * [rsvs' 1]');  % i,j,k index into vol1 for location rsvs
end

disp('---------------------------------------------------------------------------------');
disp('Initial coregistration...');

% --- Select image slices and scale for display ---
[I1,slice1,xv1,yv1,zv1] = extractslice(I1vol,r1vol,options.slicenum,options.imbright,p1(1),options.displane);
if (~options.svsmask && ~options.template_only)
    [I3,slice3,xv3,yv3,zv3] = extractslice(I3vol,r3vol,options.slicenum,options.imbright,p3(1),options.displane);
end

% --- GUI - Show the initial coregistration of the TARGET and TEMPLATE ---
if (~isempty(winid))
    axes(winid(2));  % set target plot window
    imshow(flipim(I1,p1(1),options.displane),[0 1]);
    if (~options.svsmask && ~options.template_only)
        hold on;
        contour(flipim(I3,p3(1),options.displane),contour_levels,'r');
        hold off;
    end
    text(1,5,'Coreg before','Color','w');
    drawnow;
end

disp('---------------------------------------------------------------------------------');
disp('Initial TEMPLATE and FOV/ROI...');

% --- GUI - Show overlay of Template FOV/ROI and Localizer ---
if (~isempty(winid))
    axes(winid(1));
    
    if (options.is_ROI)
        set(texthandle(1),'String',p2roi);
        [I1,slice1,xv1,yv1,zv1] = extractslice(I1vol,r1vol,ijk1roi,options.imbright,p1(1),options.displane); % For ROI, show slice instersecting middle of ROI
    elseif (is_spec)
        set(texthandle(1),'String',p2);
        [I1,slice1,xv1,yv1,zv1] = extractslice(I1vol,r1vol,ijk1svs,options.imbright,p1(1),options.displane); % For SVS, show slice instersecting middle of SVS FOV
    else
        set(texthandle(1),'String',p2);
    end
    
    % --- Paste TEMPLATE slice in 3-Space ---
    if (min(size(xv1)) < 2)
        fprintf(2,'WARNING: Cannot replane single slice Template Localizer. Skipping display.\n');
    else
        ME_warp(xv1,yv1,zv1,I1);
    end
    hold on
    
    % --- draw box around ROI found in thresholded TEMPLATE FOV ---
    if (options.is_ROI)
        plot_box(ctop2roi,cbot2roi,'b',2)
        [xr,yr,zr,xt,yt,zt] = getrange(r1vol,cbot2roi,ctop2roi,options.displane);
        
        % --- draw box showing Template FOV ---
    else
        plot_box(ctop2,cbot2,'g',2)
        plot_box(c2,c2,'g',1);  % middle slice of FOV w/ thin line
        hfill = fill3(c2(1:4,1),c2(1:4,2),c2(1:4,3),'g');
        hfill.FaceAlpha=0.2;
        [xr,yr,zr,xt,yt,zt] = getrange(r1vol,cbot2,ctop2,options.displane);
    end
    
    % --- Color in where the ROI voxels are ---
    if (options.is_ROI)
        plot3( r2vol(1,proi(:)), r2vol(2,proi(:)), r2vol(3,proi(:)), 'r.');
    end
    
    % --- Specify viewing angle ---
    setview(options.displane);
    hold off;
    axis off;
%    [xr,yr,zr,xt,yt,zt] = getrange(r1vol,options.displane);
%    axis([xr yr zr]);
    axis([xr yr zr]*1.1);
    text(xt,yt,zt,'Template Space','Color','w');
    drawnow;
    h1 = rotate3d(winid(1));    % enable rotating of view with mouse
    set(h1,'Enable','on');
end

% --- Only draw Template LOC/FOV overlay, thenn exit ---
if (options.template_only)
    fprintf(1,'Only asked to draw overlay of Template LOC & FOV. Done.\n');
    status = 1;
    return
end

% --- Generate NIFTI with SVS FOV mask on Template Loc grid ---
if (options.svsmask)
    if (options.redraw), status = 1; return; end
    
    disp('---------------------------------------------------------------------------------');
    disp('Generating SVS mask in Template Localizer space...');
    % if ((~isempty(findstr(p2,'<'))) || (~isempty(findstr(p2,'>'))))
    if ((~isempty(findstr(p2,'<'))) || (~isempty(findstr(p2,'>'))) || abs(rot2) > 1)
        fprintf(1,'Handling oblique SVS voxel...\n');
        [v1,v2,v3] = box_orient(c2);
        Vsvs_in =  [v1 v2 v3];
        Vsvs_out = [[1;0;0] [0;1;0] [0;0;1]];
        Msvs = Vsvs_out/Vsvs_in;     % this is the transform that rotates the SVS voxel into a non-oblique orientation
        cbot2x = (Msvs * cbot2')';
        ctop2x = (Msvs * ctop2')';
        r1volx = Msvs * r1vol(:,:);  % now the voxel coords of the Template Loc are rotated so the SVS voxel is non-oblique
    else
        ctop2x = ctop2;
        cbot2x = cbot2;
        r1volx = r1vol;
    end
    
    % make 3D image with mask of pixels inside SVS FOV
    Isvs = uint16(I1vol * 0);
    xr = [min([cbot2x(:,1) ;ctop2x(:,1)]) max([cbot2x(:,1) ;ctop2x(:,1)])];
    yr = [min([cbot2x(:,2) ;ctop2x(:,2)]) max([cbot2x(:,2) ;ctop2x(:,2)])];
    zr = [min([cbot2x(:,3) ;ctop2x(:,3)]) max([cbot2x(:,3) ;ctop2x(:,3)])];
    psvs =  (r1volx(1,:) >= xr(1)) & (r1volx(1,:) <= xr(2)) & ...
            (r1volx(2,:) >= yr(1)) & (r1volx(2,:) <= yr(2)) & ...
            (r1volx(3,:) >= zr(1)) & (r1volx(3,:) <= zr(2)) ;
    if (~any(psvs))
        fprintf(2,'ERROR: No voxels in Template Localizer fall inside the SVS voxel! No mask image generated.\n');
        return
    end
    Isvs(psvs) = 1;
    
    % write the Nifti file with SVS mask
    svsvol = vol1;   % copy SPM image volume of Template loc
    svsvol.fname = [svsmaskpath 'svsmask.nii'];
    svsvol.descrip = 'SVS mask';
    spm_write_vol(svsvol,Isvs);
    fprintf(1,'Saved SVS mask in %s\n',svsvol.fname);
    status = 1;
    return
end

% --- Coregister TARGET Localizer to TEMPLATE localizer ---
if (~options.reuse_xform && ~options.redraw)
    disp('---------------------------------------------------------------------------------');
    disp('Registering TARGET loc to TEMPLATE loc...');
    niftifile4  = [scratchpath 'TargetLoc2Template.nii'];
    [mat,rvol3] = nifti_coreg(vol1,vol3,options.coregmethod,niftifile4); % this returns re-sliced vol3 on vol1 grid
    fprintf(1,'Saved result in %s\n',rvol3.fname);
    last_Mat    = mat;
    last_rvol3  = rvol3;
else
    disp('---------------------------------------------------------------------------------');
    disp('Re-using previous coregistration xform betwen TARGET and TEMPLATE locs...');
    mat   = last_Mat;
    rvol3 = last_rvol3;
end

% --- Reslice the TemplateFOV/ROI into TargetLoc space (for off-line checking ONLY - not used here) ---
if (~options.redraw && ~is_spec)
    disp('---------------------------------------------------------------------------------');
    disp('Transforming Template FOV/ROI to TARGET loc space...');
    if (options.is_ROI)
        niftifile5 = [scratchpath 'TemplateROI2Target.nii'];
        rvol2 = nifti_reslice(vol3,vol2,mat,niftifile5,0); % arg5 = 0 means use nearest neighbor sampling
    else
        niftifile5 = [scratchpath 'TemplateFOV2Target.nii'];
        rvol2 = nifti_reslice(vol3,vol2,mat,niftifile5);
    end
    fprintf(1,'Saved result in %s\n',rvol2.fname);
    %    I5vol      = spm_read_vols(rvol2);
end

% --- Apply transform to TEMPLATE FOV ---
c2x    = transform_coords(c2,mat);
cbot2x = transform_coords(cbot2,mat);
ctop2x = transform_coords(ctop2,mat);

% --- Apply transform to TEMPLATE ROI bounding box ---
if (options.is_ROI)
    c2roix    = transform_coords(c2roi,mat);
    cbot2roix = transform_coords(cbot2roi,mat);
    ctop2roix = transform_coords(ctop2roi,mat);
end

% --- Move TEMPLATE ROI voxels into TARGET space ---
if (options.is_ROI)
    xyzroi  = [r2vol(:,proi(:)) ; ones(1,roicount) ];   % xyz1 coords of ROI voxels in TEMPLATE space
    xyzroix = mat * xyzroi;                             % xyz1 coords of ROI in TARGET space
    rroix   = mat * [rroi ; 1];                         % center of ROI in TARGET space
    ijk3    = round(inv(vol3.mat) * rroix);             % i,j,k index into vol3 for center of rroi (in TARGET space)
    o2roix  = mat_N2S*rroix(1:3);                       % center of mass in Siemens coords of TARGET space
end

% --- Get Siemens prescription for FOV/ROI ---
if (options.is_ROI)
    [p2x,o2x,rot2x] = box_prescrip((mat_N2S*c2roix')');
    [lx,ly,lz]      = box_size((mat_N2S*c2roix')', (mat_N2S*ctop2roix')', (mat_N2S*cbot2roix')');
else
    [p2x,o2x,rot2x] = box_prescrip((mat_N2S*c2x')');
    [lx,ly,lz]      = box_size((mat_N2S*c2x')', (mat_N2S*ctop2x')', (mat_N2S*cbot2x')');
end

disp('---------------------------------------------------------------------------------');
disp('Displaying coregistered TEMPLATE and TARGET...');

% --- Read the TargetLoc coregisted to the TemplateLoc space ---
I4vol = spm_read_vols(rvol3);

% --- Select image slices and scale for display ---
I1 = extractslice(I1vol,r1vol,options.slicenum,options.imbright,p1(1),options.displane);
I4 = extractslice(I4vol,-1,options.slicenum,options.imbright,p1(1),options.displane);

% --- GUI - Show the final coregistration of the TARGET and TEMPLATE ---
if (~isempty(winid))
    axes(winid(4));  % set target plot window
    imshow(flipim(I1,p1(1),options.displane),[0 1]);
    hold on;
    contour(flipim(I4,p1(1),options.displane),contour_levels,'r');
    hold off;
    text(1,5,'Coreg after','Color','w');
    drawnow;
end

disp('---------------------------------------------------------------------------------');
disp('TARGET and xformed and FOV/ROI...');

% --- GUI - Overlay FOV/ROI on Target Localizer ---
if (~isempty(winid))
    axes(winid(3));  % set target plot window
    set(texthandle(2),'String',p2x);
    
    if (options.is_ROI)
        [I3,slice3,xv3,yv3,zv3] = extractslice(I3vol,r3vol,ijk3,options.imbright,p3(1),options.displane);
    end
    
    % --- Paste TARGET slice in 3-Space ---
    if (min(size(xv1)) < 2)
        fprintf(2,'WARNING: Cannot replane single slice Target Localizer. Skipping display.\n');
    else
        ME_warp(xv3,yv3,zv3,I3);
    end
    hold on
    
    % --- draw box around ROI found in thresholded TEMPLATE FOV ---
    if (options.is_ROI)
        plot_box(ctop2roix,cbot2roix,'b',2)
        [xr,yr,zr,xt,yt,zt] = getrange(r3vol,cbot2roix,ctop2roix,options.displane);
        
        % --- draw box showing entire Template FOV of all slices ---
    else
        plot_box(ctop2x,cbot2x,'g',2)
        plot_box(c2x,c2x,'g',1);  % middle slice of FOV w/ thin line
        hfill = fill3(c2x(1:4,1),c2x(1:4,2),c2x(1:4,3),'g');
        hfill.FaceAlpha=0.2;
        [xr,yr,zr,xt,yt,zt] = getrange(r3vol,cbot2x,ctop2x,options.displane);
    end
    
    % --- Color in where the ROI voxels are ---
    if (options.is_ROI)
        plot3( xyzroix(1,:), xyzroix(2,:), xyzroix(3,:), 'r.');
    end
    
    % --- Specify viewing angle ---
    setview(options.displane);
    hold off;
    axis off;
%    [xr,yr,zr,xt,yt,zt] = getrange(r3vol,options.displane);
    axis([xr yr zr]);
    text(xt,yt,zt,'Target Space','Color','w');
    drawnow;
    h3 = rotate3d(winid(1));    % enable rotating of view with mouse
    set(h3,'Enable','on');
end

% --- Display the coreg matrix params ---
disp('---------------------------------------------------------------------------------');
disp('Coregistration params...');
params  = spm_imatrix(mat);
fprintf(1,'Tx,Ty,Tz = %6.2f,%6.2f,%6.2f\n',params(1),params(2),params(3));
fprintf(1,'R1,R2,R3 = %6.2f,%6.2f,%6.2f\n',rad2deg(params(4)),rad2deg(params(5)),rad2deg(params(6)));
fprintf(1,'Sx,Sy,Sz = %6.2f,%6.2f,%6.2f\n',params(7),params(8),params(9));
fprintf(1,'Ax,Ay,Az = %6.2f,%6.2f,%6.2f\n',params(10),params(11),params(12));

disp('---------------------------------------------------------------------------------');
fprintf(1,'New FOV prescription: %s (rot = %1.1f)\n',p2x,rot2x);
if (options.is_ROI), fprintf(1,'    center of mass:   L%1.1f P%1.1f H%1.1f\n',o2roix(1),o2roix(2),o2roix(3)); end
fprintf(1,'              extent: %1.1f x %1.1f x %1.1f\n',lx,ly,lz);

% --- Write new prescription to a text file, appending if exists ---
if (~isempty(dicompath))
   fp=fopen([dicompath filesep() 'imscribe_prescriptions.txt'],'a+');
   if (fp ~= -1)
       fprintf(fp,'%s   %s (rot = %1.1f)\n',datestr(now()),p2x,rot2x);
       fclose(fp);
   end
end

disp('---------------------------------------------------------------------------------');
toc;
if (options.dicom_mask && ~options.redraw)
    dcm3vol = I3vol;
    imax = max(dcm3vol(:));
    
    ijk3all = round(inv(vol3.mat) * xyzroix);       % i,j,k index into vol3 for all ROI points (in TARGET space)
    proi    = sub2ind(size(dcm3vol),ijk3all(1,:),ijk3all(2,:),ijk3all(3,:)); % convert 3 sub-indexes to 1d index (not sure why i have to do this)
    
    % --- create grayscale images w/ bright white ROI ---
    if (~options.dicom_RGB)
        roilabel = min([round(1.1*imax) 4095]);
        dcm3vol(proi) = roilabel;
        dcm3vol = flipdim(permute(dcm3vol,[2,1,3]),1);  % kluge to match data order in original dicom files (as in files3)
        dcmopts.force_UINT16 = 1;
    
    % --- create RGB images w/ RED ROI ---
    else
        dcm3vol = double(dcm3vol)/double(imax)*255;     % scale 0-255
        dcm3vol = uint8(dcm3vol);                       % RGB dicoms must be 8-bit (I think)
        redvol  = dcm3vol;
        redvol(proi)  = 255;
        dcm3vol(proi) = 0;
        dcm3vol = flipdim(permute(dcm3vol,[2,1,3]),1);  % kluge to match data order in original dicom files (as in files3)
        redvol  = flipdim(permute(redvol, [2,1,3]),1);  
        RGBvol  = cat(4,redvol,dcm3vol,dcm3vol);        % this makes NXxNYxNZx3 volume
        dcm3vol = permute(RGBvol,[1,2,4,3]);            % this makes NXxNYx3xNZ
    end    
        
    dcmopts.verbose = 1;
    dcmopts.auto_subfolder = 1;
    dcmopts.series_subfolder = 1;
%    dicom_create(files3,dcm3vol,dicompath,options.dicom_seriesnumber,['Imscribe_' options.dicom_seriesname],dcmopts);
    dicom_create(files3,dcm3vol,'./',options.dicom_seriesnumber,['Imscribe_' options.dicom_seriesname],dcmopts);
%    dicom_create(hdr3,dcm3vol,dicompath,options.dicom_seriesnumber,['Imscribe_' options.dicom_seriesname],dcmopts); % this doesn't work cuz SPM hdrs != dicominfo() hdrs

    % put a copy of the ROI in Target space NIFTI in the folder of the Target ROI
    roifilename = sprintf('Imscribe_S%4.4d_%s.nii',options.dicom_seriesnumber,options.dicom_seriesname);
    roifilepath = fileparts(files2{1});
    [stat,message,messageid] = copyfile(niftifile5,[roifilepath filesep() roifilename],'f');
    fprintf(1,'   Wrote ROI in TARGET space to %s\n',[roifilepath filesep() roifilename]);    
end
disp('---------------------------------------------------------------------------------');
disp('Done.');

status = 1;
return;


% -------------------------------------------------------------------------
% SUBROUTINES
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function plot_box(ctop,cbot,color,width)
if (nargin < 3), color = 'k'; end
if (nargin < 4), widhth = 1; end

plot3(cbot(:,1),cbot(:,2),cbot(:,3), color,'LineWidth',width);
plot3(ctop(:,1),ctop(:,2),ctop(:,3), color,'LineWidth',width);
for j=1:4, plot3([cbot(j,1);ctop(j,1)], [cbot(j,2);ctop(j,2)], [cbot(j,3);ctop(j,3)],color,'LineWidth',width); end
return

% -------------------------------------------------------------------------
% Flip a 2D image so it displays correctly w/ imshow
function [imout] = flipim(imin,orient_code,plane)

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
return

% -------------------------------------------------------------------------
% Pull slice from 3D volume. Get desired plane given orientation of 3D input vol
function [im,s0,xv,yv,zv] = extractslice(vol,r,slice,bright,orient_code,plane)
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
return

% -------------------------------------------------------------------------
% clip input array to min,max limits
%  (uses matrix Logicals to replace Loops)
function  y = clip( x, lo, hi )
y = (x .* [x<=hi])  +  (hi .* [x>hi]);
y = (y .* [x>=lo])  +  (lo .* [x<lo]);
return

% -------------------------------------------------------------------------
% Get row,column,normal length of FOV bounding box
function [lx,ly,lz] = box_size(c,ctop,cbot)
lx = sqrt( (c(2,1)-c(1,1))^2 + (c(2,2)-c(1,2))^2 + (c(2,3)-c(1,3))^2 );
ly = sqrt( (c(3,1)-c(2,1))^2 + (c(3,2)-c(2,2))^2 + (c(3,3)-c(2,3))^2 );
lz = sqrt( (ctop(1,1)-cbot(1,1))^2 + (ctop(1,2)-cbot(1,2))^2 + (ctop(1,3)-cbot(1,3))^2 );
return

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
return

% -------------------------------------------------------------------------
function [xr,yr,zr,xt,yt,zt] = getrange(rvol,cbot,ctop,plane)
xr = [min(rvol(1,:)) max(rvol(1,:))];  % Xrange extent of vol
yr = [min(rvol(2,:)) max(rvol(2,:))];
zr = [min(rvol(3,:)) max(rvol(3,:))];

xr = [min([cbot(:,1) ;ctop(:,1); xr']) max([cbot(:,1) ;ctop(:,1); xr'])]; % Xrange extent of box and vol
yr = [min([cbot(:,2) ;ctop(:,2); yr']) max([cbot(:,2) ;ctop(:,2); yr'])];
zr = [min([cbot(:,3) ;ctop(:,3); zr']) max([cbot(:,3) ;ctop(:,3); zr'])];

if (xr(1) == xr(2)), xr(1) = xr(1)-.5; xr(2) = xr(2)+.5; end % kluge for single slice so axis(xr,yr,zr) doesn't crash
if (yr(1) == yr(2)), yr(1) = yr(1)-.5; yr(2) = yr(2)+.5; end
if (zr(1) == zr(2)), zr(1) = zr(1)-.5; zr(2) = zr(2)+.5; end

switch plane % return location to insert text on top of pictures
    case 1      % Axial slice
        xt = xr(2); yt = yr(2)-10; zt = zr(2)-10;
    case 2      % Coronal slice
        xt = xr(2); yt = yr(1); zt = zr(2)-10;
    case 3      % Sagittal slice
        xt = xr(2); yt = yr(2); zt = zr(2)-10;
end
return



