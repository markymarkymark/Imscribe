% to compile
%
%  mcc -m imscribe -a spm_dicom_dict.mat 

function varargout = imscribe(varargin)
% IMSCRIBE M-file for imscribe.fig
%      IMSCRIBE, by itself, creates a new IMSCRIBE or raises the existing
%      singleton*.
%
%      H = IMSCRIBE returns the handle to a new IMSCRIBE or the handle to
%      the existing singleton*.
%
%      IMSCRIBE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMSCRIBE.M with the given input arguments.
%
%      IMSCRIBE('Property','Value',...) creates a new IMSCRIBE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imscribe_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imscribe_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imscribe

% Last Modified by GUIDE v2.5 11-Jul-2019 12:13:54

% --- initialize imscribe environment ---
%stat = imscribe_init(); % moved this to imscribe_OpeningFcn()
%if (~stat), return; end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imscribe_OpeningFcn, ...
                   'gui_OutputFcn',  @imscribe_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before imscribe is made visible.
function imscribe_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imscribe (see VARARGIN)

% Choose default command line output for imscribe
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imscribe wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- initialize GUI objects ---
imscribe_reset_gui(handles);

% --- initialize imscribe environment ---
stat = imscribe_init(); % moved this call to main routine
%if (~stat), buttonEXIT_Callback(hObject, eventdata, handles); end    % this Fails - wrong place for exiting GUI
if (~stat), error('terminating program'); end    



% --- Outputs from this function are returned to the command line.
function varargout = imscribe_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in textPRESCRIP_TARGET.
function textPRESCRIP_TARGET_Callback(hObject, eventdata, handles)
% hObject    handle to textPRESCRIP_TARGET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of textPRESCRIP_TARGET





% --- Executes on button press in buttonTEMPLATE.
function buttonTEMPLATE_Callback(hObject, eventdata, handles)
% hObject    handle to buttonTEMPLATE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global path1 path1files options

% --- get user specified file/files ---
[path1,path1files,is_nifti,nfiles,is_rda] = imscribe_choose_files(1);
if (isequal(path1,'')), return; end
if (is_rda)
    fprintf(2,'ERROR - RDA files should only be selected for the TEMPLATE FOV/ROI\n');
    path1 = '';
    return
end
options.is_nifti(1) = is_nifti;

% --- display folder for multiple files, filename for Nifti ---
if (nfiles == 1), dispname = path1files{1};
else              dispname = path1; end
slen = 60;  % number of characters that fit in text box
l    = length(dispname);
if (l > 0)
    if (l > slen), dispname = ['...' dispname(l-slen:l)]; end % can only fit so many chars in widget
    set(handles.textTEMPLATE_FOLDER,'string',dispname);
end

% --- Can't re-use previous coreg xform if new TEMPLATE or TARGET ---
set(handles.checkbox_REUSE_xform,'value',0);

return

% --- Executes on button press in buttonTEMPLATE_ROI.
function buttonTEMPLATE_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to buttonTEMPLATE_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global path2 path2files options

% --- get user specified file/files ---
[path2,path2files,is_nifti,nfiles,is_rda] = imscribe_choose_files(2);
if (isequal(path2,'')), return; end
options.is_nifti(2) = is_nifti;
options.is_rda      = is_rda;

% --- display folder for multiple files, filename for Nifti ---
if (nfiles == 1), dispname = path2files{1};
else              dispname = path2; end
slen = 60;  % number of characters that fit in text box
l    = length(dispname);
if (l > 0)
    if (l > slen), dispname = ['...' dispname(l-slen:l)]; end % can only fit so many chars in widget
    set(handles.textTEMPLATE_ROI_FOLDER,'string',dispname);
end

return

% --- Executes on button press in buttonTARGET.
function buttonTARGET_Callback(hObject, eventdata, handles)
% hObject    handle to buttonTARGET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global path3 path3files options

% --- get user specified file/files ---
[path3,path3files,is_nifti,nfiles,is_rda] = imscribe_choose_files(3);
if (isequal(path3,'')), return; end
if (is_rda)
    fprintf(2,'ERROR - RDA files should only be selected for the TEMPLATE FOV/ROI\n');
    path3 = '';
    return
end
options.is_nifti(3) = is_nifti;

% --- display folder for multiple files, filename for Nifti ---
if (nfiles == 1), dispname = path3files{1};
else              dispname = path3; end
slen = 60;  % number of characters that fit in text box
l    = length(dispname);
if (l > 0)
    if (l > slen), dispname = ['...' dispname(l-slen:l)]; end % can only fit so many chars in widget
    set(handles.textTARGET_FOLDER,'string',dispname);
end

% --- Can't re-use previous coreg xform if new TEMPLATE or TARGET ---
set(handles.checkbox_REUSE_xform,'value',0);
return

% --- Executes on button press in buttonGET_RAW.
function buttonGET_RAW_Callback(hObject, eventdata, handles)
% hObject    handle to buttonGET_RAW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

slen = 60;  % number of characters that fit in text box

path = imscribe_choose_raw();
l    = length(path);
if (l > 0)
    if (l > slen), path = ['...' path(l-slen:l)]; end % can only fit so many chars in widget
    set(handles.textTARGET_FOLDER,'string',path);
end

% --- Executes on button press in buttonGO.
function buttonGO_Callback(hObject, eventdata, handles, do_redraw )
% hObject    handle to buttonGO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global scratchpath dicom_path
global path1 path2 path3 path1files path2files path3files
global options

%  which button was pushed?
if (nargin < 4), do_redraw  = 0; end
options.redraw  = do_redraw;
options.svsmask = get(handles.buttonSVSMASK,'Value');
options.hippo_slab  = get(handles.radiobuttonHIPPOSLAB,'Value');
options.template_only = get(handles.buttonTEMPLATES_ONLY,'Value');

% make sure paths are specified
if (options.svsmask)
    if (isempty(path1) || isempty(path2))
        errordlg('ImScribe: For SVSMASK, you must specify the TEMPLATE files');
        return;
    end
elseif (options.template_only)
    if (isempty(path1) || isempty(path2))
        errordlg('ImScribe: For DRAW TEMPLATES, you must specify both TEMPLATE files');
        return;
    end
elseif (options.hippo_slab)
    if (isempty(path3))
        errordlg('ImScribe: For Hippo Slab, you must specify the TARGET Localizer');
        return;
    end
else
    if (isempty(path1) || isempty(path2) || isempty(path3))
        errordlg('ImScribe: You must specify the TEMPLATE and TARGET files');
        return;
    end
end   

% setup for call to imscribe_engine() 
paths{1}   = scratchpath; % scratch folder
paths{2}   = scratchpath; % SVS mask result folder
paths{3}   = dicom_path;  % Dicom mask result folder
winid      = [handles.window1 handles.window2 handles.window3 handles.window4]; 
texthandle = [handles.textPRESCRIP_TEMPLATE handles.textPRESCRIP_TARGET];
if     (get(handles.radioRIGID,     'Value'))
    options.coregmethod = 'rigidbody';
elseif (get(handles.radioAFFINE,    'Value'))
    options.coregmethod = 'affine';
elseif (get(handles.radioINTERMODAL,'Value'))
    options.coregmethod = 'coreg';
end
options.displane     =  1 * isequal(get(handles.radioTRANSVERSE,'Value'),1) ...
                      + 2 * isequal(get(handles.radioCORONAL,   'Value'),1) ...
                      + 3 * isequal(get(handles.radioSAGITTAL,  'Value'),1);
options.slicematch  = 0;
options.contourlevel= get(handles.sliderCONTOUR,'Value');
options.slicenum    = get(handles.sliderSLICENUM,'Value');
options.imbright    = get(handles.sliderBRIGHTNESS,'Value');
options.is_ROI      = get(handles.checkboxROI,'Value');
options.ROIthresh   = str2num(get(handles.textinput_ROIthresh,'String'));
options.reuse_xform = get(handles.checkbox_REUSE_xform,'Value');
options.dicom_mask  = get(handles.buttonDICOMMASK,'Value');
options.dicom_RGB   = 0;
options.skullstrip_temploc = get(handles.checkbox_SSTempLoc,'Value');
options.skullstrip_targloc = get(handles.checkbox_SSTargLoc,'Value');

% make the button that was pushed be RED 
bcolor = get(handles.buttonGO,'BackGroundColor');
if     (options.redraw),  set(handles.pushbuttonREDRAW,'BackGroundColor',[1 0 0]);
else                      set(handles.buttonGO,        'BackGroundColor',[1 0 0]); end
watchon;
if (~options.hippo_slab)    
    imscribe_engine(paths,path1files,path2files,path3files,options,winid,texthandle);
else
    imscribe_custom(paths,path1files,path2files,path3files,options,winid,texthandle); % Handle custom actions
end
    
watchoff;
set(handles.pushbuttonREDRAW,'BackGroundColor',bcolor);
set(handles.buttonGO,'BackGroundColor',bcolor);
set(handles.buttonDICOMMASK,'Value',0);  % always turn off dicom write option after completion
return

% --- Executes during object deletion, before destroying properties.
function sliderSLICENUM_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to sliderSLICENUM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in buttonEXIT.
function buttonEXIT_Callback(hObject, eventdata, handles)
% hObject    handle to buttonEXIT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1);




% --- Executes on button press in checkboxROI.
function checkboxROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxROI



function textinput_ROIthresh_Callback(hObject, eventdata, handles)
% hObject    handle to textinput_ROIthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textinput_ROIthresh as text
%        str2double(get(hObject,'String')) returns contents of textinput_ROIthresh as a double


% --- Executes during object creation, after setting all properties.
function textinput_ROIthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textinput_ROIthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_RESETviews.
function pushbutton_RESETviews_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RESETviews (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- reset the rotation views of the image windows ---
plane =  1 * isequal(get(handles.radioTRANSVERSE,'Value'),1) ...
       + 2 * isequal(get(handles.radioCORONAL,   'Value'),1) ...
       + 3 * isequal(get(handles.radioSAGITTAL,  'Value'),1);
switch plane
    case 1      % Axial slice
        angles=[180,90];
    case 2      % Coronal slice
        angles=[180,0];
    case 3      % Sagittal slice
        angles=[89.9,0.1];
end
view(handles.window1,angles);
view(handles.window3,angles);
		


% --- Executes on button press in checkbox_REUSE_xform.
function checkbox_REUSE_xform_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_REUSE_xform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_REUSE_xform


% -----------------------SLIDERS------------------------------------

% --- Executes on slider movement.
function sliderBRIGHTNESS_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSLICENUM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderBRIGHTNESS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSLICENUM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function sliderSLICENUM_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSLICENUM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderSLICENUM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSLICENUM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function sliderCONTOUR_Callback(hObject, eventdata, handles)
% hObject    handle to sliderCONTOUR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderCONTOUR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderCONTOUR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% -----------------------RADIO BUTTONS------------------------------------

% --- Executes on button press in radioAFFINE.
function radioAFFINE_Callback(hObject, eventdata, handles)
% hObject    handle to radioAFFINE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioAFFINE


% --- Executes on button press in radioRIGID.
function radioRIGID_Callback(hObject, eventdata, handles)
% hObject    handle to radioRIGID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioRIGID

% --- Executes on button press in buttonSVSMASK.
function buttonSVSMASK_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSVSMASK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Do the same thing as the GO button, but set option to generate mask of SVS voxel on Template loc ---
%buttonGO_Callback(hObject, eventdata, handles, 0,1);


% --- Executes on button press in buttonDICOMMASK.
function buttonDICOMMASK_Callback(hObject, eventdata, handles)
% hObject    handle to buttonDICOMMASK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttonDICOMMASK
global options
persistent series_num

if (isempty(series_num)), series_num = 500; end

if (get(hObject,'Value'))
    prompt = {'Series name:','Series number:'};
    dlg_title = 'Dicom Mask params';
    num_lines = 1;
    def = {'roi1',sprintf('%1d',series_num)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if (isempty(answer))
        set(hObject,'Value',0); 
    else
        options.dicom_seriesname = answer{1};
        options.dicom_seriesnumber = str2num(answer{2});
        series_num = options.dicom_seriesnumber + 1;
    end
end

% --- Executes on button press in radiobuttonHIPPOSLAB.
function radiobuttonHIPPOSLAB_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonHIPPOSLAB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonHIPPOSLAB
return

% --- Executes on button press in radiobuttonHIPPOSLAB.
function buttonTEMPLATES_ONLY_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonHIPPOSLAB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonHIPPOSLAB
return

% --- Executes on button press in pushbuttonREDRAW.
function pushbuttonREDRAW_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonREDRAW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Do the same thing as the GO button, but set option to just redraw the pics ---
buttonGO_Callback(hObject, eventdata, handles, 1);
return


% --- Executes on button press in buttonAUTOLOC.
function buttonAUTOLOC_Callback(hObject, eventdata, handles)
% hObject    handle to buttonAUTOLOC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global path3 path3files RT_export_path options

%msgbox('Searching for most recent RTexport dicom series.\nPlease wait...','ImScribe','modal'); 

% Get most recent RTexport folder
dirs       = dir(RT_export_path);
dirs       = dirs(3:numel(dirs));   % skip '.' and '..'
dirindex   = [dirs.isdir];
dirs       = dirs(dirindex == 1);   % weed out non directories
dates      = [dirs.datenum];
[tmp,i]    = sort(-dates);
newest_dir = dirs(i(1)).name;
path3      = [RT_export_path filesep() newest_dir filesep()];
files      = dir([path3 '*.dcm']);
nfiles     = numel(files);
if (nfiles == 0)
    msg = sprintf('Most recent RTexport folder:\n\n%s\n\nContains no Dicoms files!',path3);
    warndlg(msg,'CreateMode','modal'); 
    uiwait(gcf);    
    path3 = '';
    set(handles.textTARGET_FOLDER,'string','');
return
end

% find all dicoms from latest series
nums   = sscanf(files(nfiles).name,'001_%d_%d.dcm'); % get series number from last filename
serstr = sprintf('001_%06d_*.dcm',nums(1));          % filter for series number
files  = dir([path3 serstr]);
nfiles = numel(files);
path3files = [];
for i=1:nfiles
    path3files{i} = [path3 files(i).name];
end
options.is_nifti(3) = 0;
msg = sprintf('Loaded %1d Dicom files from series %1d in \n\n%s',nfiles,nums(1),path3);
msgbox(msg,'ImScribe','modal'); 
uiwait(gcf);   
    
dispname = path3;
slen = 60;  % number of characters that fit in text box
l    = length(dispname);
if (l > 0)
    if (l > slen), dispname = ['...' dispname(l-slen:l)]; end % can only fit so many chars in widget
    set(handles.textTARGET_FOLDER,'string',dispname);
end

% --- Can't re-use previous coreg xform if new TEMPLATE or TARGET ---
set(handles.checkbox_REUSE_xform,'value',0);

return

% --- Executes on button press in checkbox_SSTempLoc.
function checkbox_SSTempLoc_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SSTempLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SSTempLoc


% --- Executes on button press in checkbox_SSTargLoc.
function checkbox_SSTargLoc_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SSTargLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SSTargLoc

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% User-written routines that interect with the GUI
% (above here is autogenerated code from the GUIDE app)
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


function imscribe_reset_gui(handles)

set(handles.textPRESCRIP_TEMPLATE,'string','');
set(handles.textPRESCRIP_TARGET,  'string','');
set(handles.textTEMPLATE_FOLDER,  'string','');
set(handles.textTEMPLATE_ROI_FOLDER,  'string','');
set(handles.textTARGET_FOLDER,    'string','');
set(handles.radioAFFINE,          'Value',1);
set(handles.radioSAGITTAL,        'Value',1);
set(handles.radiobuttonSTANDARD,  'Value',1);

axis(handles.window1,'off');
axis(handles.window2,'off');
axis(handles.window3,'off');
axis(handles.window4,'off');

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% User-written routines that DO NOT interect with the GUI
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


function status = imscribe_init()

global RT_export_path scratchpath dicom_path
global path1 path2 path3
global options
global last_Mat

status = 0;

% initialize Nifti/Dicom flag for the 2 image sets
options.is_nifti = [0 0 0];
options.is_rda   = 0;   % in place of SVS for voxel mask

% discard any record of last xformation matrix found by imscribe_engine()
last_Mat = [];

% reset paths for input images
path1 = '';
path2 = '';
path3 = '';

% set the default folders to look in 
[p1,p2,p3,p4,p5,p6] = imscribe_default_folders();
if (isempty(p4)), return; end
scratchpath            = p4;
RT_export_path         = p5; 
dicom_path             = p6; 

% Initialize default folders for file selection ---
imscribe_choose_files({p1,p2,p3});

status = 1;
return

% -------------------------------------------------------------------------
function [tpath,file,is_nifti,nfiles,is_rda] = imscribe_choose_files(slot)
% Prompt user to select the LOCALIZER or TARGET files 
% Can be multiple Dicoms or single Dicom or single Nifti

persistent last_paths       % these are the previously user-selected file paths 

tpath    = '';
file     = '';
is_nifti = 0;
is_rda   = 0;
nfiles   = 0;

% --- initialize default file path path for each use (aka "slot") ---
if (iscell(slot))
    last_paths = slot;
    return
end

% --- make sure this routine was initialized ---
if (isempty(last_paths)), error('imscribe_choose_files() routine was not initialized!'); end
    
% --- Select all needed files  ---
%wildcard   = {'*.dcm;*.IMA','Dicoms (*.dcm, *.IMA)'; '*.nii;*.nii.gz','Niftis (*.nii, *.nii.gz)'};
%[file, tpath] = uigetfile(wildcard,'Select one NIFTI or ALL desired Dicom files',last_paths{slot},'MultiSelect','on');
file = me_uipickfiles('FilterSpec',last_paths{slot},'Prompt','Select one NIFTI or ALL desired Dicom files','REFilter','\.dcm$|\.nii$|\.IMA$\.MR$\.nii.gz$|\.rda$');
% 'Type',{'*.dcm','Dicoms';'*.nii','NIFTIs'}
if isequal(file,0), return; end % user chose cancel
if isempty(file),   return; end % user chose done w/ no files selected!!

% --- get path for next time using this UI ---
tpath = fileparts(file{1});
last_paths{slot} = tpath; % remember for next time

nfiles = numel(file);
if nfiles > 1, file = sort(file); end   % sort filenames in ascending order in case converting dicoms needs it!

% --- Figure out if Dicom or NIFTI or .rda ---
[~, ~, ext] = fileparts(file{1});
if (isequal(ext,'.dcm'))
    is_nifti = 0;
elseif (isequal(ext,'.IMA'))
    is_nifti = 0;
elseif (isequal(ext,'.MR'))
    is_nifti = 0;
elseif (isequal(ext,'.nii') || isequal(ext,'.gz'))
    if (nfiles > 1)
        tpath = '';
        fprintf(2,'ERROR - Only one NIFTI file should be selected\n',ext);
    else
        is_nifti = 1;
    end
elseif (isequal(ext,'.rda'))
    if (nfiles > 1)
        tpath = '';
        fprintf(2,'ERROR - Only one RDA file should be selected\n',ext);
    else
        is_rda = 1; 
    end
else
    tpath = '';
    fprintf(2,'ERROR - unrecognized file name extension: %s\n',ext);
end
if (isequal(ext,'.gz'))
        fprintf(1,'Unzipping %s...\n',file{1});
        file = gunzip(file); % need to uncompress .gz files for SPM to read
end
return


