function [fid,hdr] = readRDA2(infile)

% --- Remember path of last time this routine was called ---
global RDA_defpath
if (size(RDA_defpath,1) == 0) 
	RDA_defpath = pwd();
end

fid = [];
hdr = [];

if (nargin < 1 || isempty(infile))
    [file1, path1] = uigetfile('*.rda','Select a Siemens RDA file',RDA_defpath);
	if isequal(file1,0)
        hdr = [];
		disp('Program cancelled.');
		return
	end
	RDA_defpath = path1;
    infile = [path1 filesep() file1];
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
fid = reshape(fid, 2, hdr.VectorSize, hdr.CSIMatrixSize0, hdr.CSIMatrixSize1, hdr.CSIMatrixSize2);
fid = complex(fid(1,:,:,:,:),fid(2,:,:,:,:));
%fid = squeeze(fid);  % this doesn't remove leading singleton!!
fidx = fid(1,:,:,:,:);
fclose(fp);

end