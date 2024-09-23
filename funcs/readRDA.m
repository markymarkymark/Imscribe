function [fid,hdr] = readRDA(infile)
% readRDA() - Read spectroscopy RDA file
%           - Return all the header fields in a struct
%
% Syntax:
%   [fid,hdr] = readRDA([rdafile])
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

fid = [];
hdr = [];
if (nargin < 1 || isempty(infile))
    infile = uigetfile_plus(mfilename(),'rdapath',pwd(),'*.rda','Select a Siemens RDA file');
	if isempty(infile); return;	end
end
fp = fopen(infile,'r');

head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';
tline = fgets(fp);
if (isempty(strfind (tline , head_start_text)))
    fprintf(2,'ERROR: First line of RDA file should contain "%s"\n',head_start_text);
    fclose(fp);
    return
end

% read header
tline = fgetl(fp);                              % Note: fgetl() removes newline at end
while (isempty(strfind(tline , head_end_text)))
    parts = regexp(tline,':','split','once');   % get name and value separated by ':'
    parts{1} = strrep(parts{1},'[', '');        % remove [ 
    parts{1} = strrep(parts{1},']', '');
    num = str2num(parts{2});                    % try to convert to a number
    if (isempty(num))                                                   % STRING variable
        cmd = sprintf('hdr.%s = ''%s'';',parts{1},strtrim(parts{2}));    
    else                                                                % NUMERIC variable 
        cmd = sprintf('hdr.%s = %s;',parts{1},parts{2});     
        cmd = erase(cmd,',');                                           % this handles VARIABLE[x,y] = ... (makes VARIABLExy = ...)
    end
%    disp(cmd);
    eval(cmd)
    tline = fgetl(fp);
end

% read data
fid = fread(fp , hdr.CSIMatrixSize0 * hdr.CSIMatrixSize1 * hdr.CSIMatrixSize2 * hdr.VectorSize * 2 , 'double');  
fclose(fp);
fid = reshape(fid, 2, hdr.VectorSize, hdr.CSIMatrixSize0, hdr.CSIMatrixSize1, hdr.CSIMatrixSize2);
fid = complex(fid(1,:,:,:,:),fid(2,:,:,:,:));
%fid = squeeze(fid);  % this doesn't remove leading singleton!!
%fidx = fid(1,:,:,:,:);
end