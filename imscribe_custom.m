% imscribe_custom() - handle for executing custom actions
function [status] = imscribe_custom(paths,files1,files2,files3,options,winid,texthandle)

% --- Set return variables ---
status = 0;	% this means error

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

% --- Start counting program duration ---
tic;

%--- TARGET LOCALIZER ---
niftifile3 = [scratchpath 'TargetLoc.nii'];
if (~options.svsmask)
    if (~options.is_nifti(3))
        fprintf(1,'   Converting %1d TARGET Localizer Dicoms to Nifti\n',size(files3,iscell(files3)+1));
        %            [vol3,hdr3] = dicom_spm2nifti(files3,niftifile3);
        vol3 = dicom_spm2nifti(files3,niftifile3);
    else
        fprintf(1,'   Reading TARGET Localizer Nifti header\n',numel(files1));
        vol3 = spm_vol(files3);
        [stat,message,messageid] = copyfile(files3,niftifile3,'f');
        fprintf(1,'   Copied TARGET Localizer to %s\n',niftifile3);
    end
end

% --- Run custom script ---
cmd = ['hippo_slab.sh ' niftifile3];
fprintf(1,'   Executing command: %s\n',cmd);
[stat, result] = system(cmd);
if (stat ~= 0)
    fprintf(2,'   %s\n',result);
    return
end
nums = sscanf(result,'%f',6);
origin = nums(1:3);
v3     = nums(4:6);

% --- Get Siemens slice/slab prescription ---
center = origin_prescrip(mat_N2S * origin);
plane  = vector_prescrip(mat_N2S * v3);
prescrip = [plane ' ' center];

% --- Write new prescription to a text file, appending if exists ---
disp('---------------------------------------------------------------------------------');
fprintf(1,'New FOV prescription: %s\n',prescrip);
if (~isempty(dicompath))
   fp=fopen([dicompath filesep() 'imscribe_prescriptions.txt'],'a+');
	if (fp > 0)
   		fprintf(fp,'%s   %s \n',datestr(now()),prescrip);
   		fclose(fp);
	else
		fprintf(2,'Cannot write new prescription to Scanner Host. You will need to copy it from above\n');
	end
end

status = 1;
return



