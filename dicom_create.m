function filenames = dicom_create(dcmfiles,data,outfolder,new_series,new_description,options,xtrahdr)
% Create new dicom file(s) for a 2D or 3D data set.
% Requires template dicom filename(s), or header structure(s) as returned by dicominfo().
%
% INPUTS:
%   dcmfiles    - string containing a valid dicom filename, or a dicom header structure (for 2D data set)
%               - a cell array of dicom filenames, or dicom header structures (for 3D data set, i.e. single series)
%   data        - a 2D or 3D image data array to put in the new dicom file(s)
%   outfolder   - the folder location to hold the new dicom file(s)
%   new_series  - the series number to assign to the new dicom files (should be differeny than the original)
%   new_description - character string with the new dicom series description
%   options     - a structure of optional features (see Example below)
%   xtrhdr      - a structure with additional dicom header fields to be changed in the output dicoms
%                 e.g. xtra.PatientName.GivenName = 'John';
%                      xtra.PatientName.FamilyName = 'Doe';
%
% OUTPUTS:
%   filenames   - cell array of newly created dicom files
%
% NOTES:
%   1) template dicom files or headers must be from a single series (for 3D data)
%   2) data set should be scaled to UINT16 for best compatiblilty
%   3) data set should be oriented the same as the original dicom data (watch out for tranposing)
%   4) To write a 3D volume to a dicom series, pass in a cell array of dicom filenames or header structs
%
% Example:
%   h = dicominfo('sample.dcm');
%   im = dicomread(h);
%   im = max(im(:))-im;     % invert image
%   options.verbose = 1;
%   options.force_UINT16 = 1;
%   options.auto_subfolder = 1;
%   newfilename = dicom_create(h,im,pwd(),500,'Inverted_image',options);
%       or
%   newfilename = dicom_create('sample.dcm',im,pwd(),500,'Inverted_image',options);
%
% MODIFICATION HISTORY:
%   Written by M.Elliott 10/2013

% --- Returned params --
filenames = {};

% --- optional params ---
if (nargin < 5)
    error('%s\n%s\n','Wrong number of input arguments.','Usage:  filenames = dicom_create(dcmfiles,data,outfolder,new_series,new_description,[options],[xtrahdr])');
end
if (nargin < 6), options.dummy = 0; end % forces options to be a structure
if (nargin < 7), xtrahdr = []; end
if (~isfield(options,'force_UINT16')),     options.force_UINT16     = 0; end
if (~isfield(options,'force_UINT8')),      options.force_UINT8      = 0; end
if (~isfield(options,'verbose')),          options.verbose          = 0; end
if (~isfield(options,'auto_subfolder')),   options.auto_subfolder   = 0; end
if (~isfield(options,'series_subfolder')), options.series_subfolder = 0; end
if (~isfield(options,'MFsort_naming')),    options.MFsort_naming    = 0; end


% --- Parse optional header fields to set in output dicoms ---
if (~isempty(xtrahdr))
    xtratags = fieldnames(xtrahdr);
    nxtra    = numel(xtratags);
else
    nxtra    = 0;
end

% --- cell array passed in - find out if filenames or dicominfo() structures ---
if (iscell(dcmfiles))
    nfiles = numel(dcmfiles);
    if (isstruct(dcmfiles{1})) headers_in = 1;         % passed in dicominfo() headers, not filenames
    else                       headers_in = 0; end     % passed in dicom filenames
    
    % --- single filename or header struct passed in ---
else
    nfiles = 1;
    if (isstruct(dcmfiles)), headers_in = 1;         % passed in a dicominfo() header, not a filename
    else                     headers_in = 0; end     % passed in a dicom filename
end

% --- Set output file location ---
if (isempty(outfolder)), outfolder = '.'; end
if (exist(outfolder,'dir') ~= 7)
    error('%s is not a valid folder name.\n',outfolder);
end

% --- Get stuff needed for output dicoms ---
if (options.force_UINT16), data = uint16(data); end
if (options.force_UINT8),  data = uint8(data); end
dmax = max(data(:));
uid  = dicomuid();    %  new dicom series needs a new Dicom UID

% --- Create new dicom files ---
for i=1:nfiles
    
    % --- get next dicom header struct ---
    if (headers_in)
        if (nfiles == 1), h = dcmfiles;
        else              h = dcmfiles{i}; end
        
        % --- read next template dicom header ---
    else
        if (nfiles == 1), h = dicominfo(dcmfiles);
        else              h = dicominfo(dcmfiles{i}); end
    end
    
    if (i == 1)
        % --- Get input image data size ---
        is_RGB  = 0;
        nd = ndims(data);
        [dx,dy,dp,dt] = size(data);
        switch (nd)
            case 2
                dz = 1;                     % this is a single 2D image
            case 3
                if (nfiles == 1)
                    if (dp ~= 1 && dp ~= 3)
                        error('Input image data is %1dx%1dx%1d. For a single dicom, image data must be 2D or 2Dx3 (i.e. RGB)\n',dx,dy,dp);
                    else
                        is_RGB = (dp == 3); % RGB data if NXxNYx3
                        dz     = 1;         % one 2D or RGB image
                    end
                else
                    dz = dp;                % this is multiple 2D images
                end
            case 4
                if (dp ~= 3)
                    error('Input image data is %1dx%1dx%1d%1d. Expected this to be multiple RGB images (i.e %1dx%1dx3xNZ)\n',dx,dy,dp,dt,dx,dy);
                else
                    is_RGB = 1;             % this is multiple RGB images
                    dz     = dt;
                end
            otherwise
                error('Input image data has %1d dimensions.\n',nd);
        end
        
        % --- Check that data size matches template dicom image size ---
        nx = h.Rows; ny = h.Columns; nz = nfiles;
        if ((dx ~= nx) || (dy ~= ny) || (dz ~= nz))
            error('Input data array size (%1d, %1d, %1d) does not match dicom image data size (%1d, %1d, %1d)\n',dx,dy,dz,nx,ny,nz);
        end
        
        % --- Sub folders ---
        if (options.auto_subfolder)
            patname = [clean_string(h.PatientName.FamilyName) '_' clean_string(checkstruct(h.PatientName,'GivenName','xxx'))];
            outfolder = [outfolder filesep() patname '_' h.StudyDate];
            if (~isdir(outfolder)), mkdir(outfolder); end
        end
        if (options.series_subfolder)
            outfolder = [outfolder filesep() sprintf('S%4.4d_%s',new_series,new_description)];
            if (~isdir(outfolder)), mkdir(outfolder); end
        end
        if (options.verbose), fprintf(1,'Creating %1d Dicom files from series %1d in %s\n',nfiles,new_series,outfolder); end
        if (options.verbose && is_RGB), fprintf(1,'Dicoms will be RGB!\n'); end
    end
    
    % --- Fill in mandatory new header values ---
    h.SeriesInstanceUID = uid;
    h.SeriesNumber      = new_series;
    h.SeriesDescription = new_description;
    h.WindowCenter      = dmax/2;
    h.WindowWidth       = dmax;
    
    % --- Fill in optional additional header items ---
    for j=1:nxtra
        h = setfield(h,xtratags{j},getfield(xtrahdr,xtratags{j}));
    end
    
    % --- Set output filename ---
    if (options.MFsort_naming)
        filenames{i} = sprintf('%s%sE1S%4.4dX01I%5.5d.dcm',outfolder,filesep(),h.SeriesNumber,h.InstanceNumber);
    else
        filenames{i} = sprintf('%s%sE1_%6.6d_%6.6d.dcm',outfolder,filesep(),h.SeriesNumber,h.InstanceNumber);
    end
    
    % --- Write the new dicom ---
    if (options.verbose), fprintf(1,'.'); end
    %    dicomwrite(data(:,:,i), filenames{i}, h, 'CreateMode', 'copy');
    if (~is_RGB)
%        dicomwrite(data(:,:,i), filenames{i}, h, 'CreateMode', 'copy', 'WritePrivate', 1, 'Dictionary', 'ME_dicom-dict.txt'); % need custom dict b/c of bug with some Matlab versions and 0B pixeldata
        dicomwrite(data(:,:,i), filenames{i}, h, 'CreateMode', 'copy', 'WritePrivate', 1);
    else
        if (~isequal(class(data),'uint8'))
            fprintf(2,'WARNING: Writing RGB data that is NOT UINT8. Might not work. Try options.force_UINT8=1.\n');
        end
        h.PhotometricInterpretation = 'RGB';
        h.SamplesPerPixel = 3;
        dicomwrite(data(:,:,:,i), filenames{i}, h, 'CreateMode', 'copy', 'WritePrivate', 1, 'Dictionary', 'ME_dicom-dict.txt'); 
    end
end
if (options.verbose), fprintf(1,'\n'); end

return

