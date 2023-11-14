function [p1,p2,p3,p4,p5,p6] = imscribe_default_folders()
% set default folders for Imscribe to find image files
% SHOULD be edited for installation-specific folder names

p1 = ''; p2 = ''; p3 = ''; p4 = ''; p5 = ''; p6 = ''; 

% --- Check saved MATLAB preferences ---
if (ispref(mfilename(),'IMSCRIBE_P1')),  p1 = getpref(mfilename(),'IMSCRIBE_P1');     end
if (ispref(mfilename(),'IMSCRIBE_P2')),  p2 = getpref(mfilename(),'IMSCRIBE_P2');     end
if (ispref(mfilename(),'IMSCRIBE_P3')),  p3 = getpref(mfilename(),'IMSCRIBE_P3');     end
if (ispref(mfilename(),'IMSCRIBE_P4')),  p4 = getpref(mfilename(),'IMSCRIBE_P4');     end
if (ispref(mfilename(),'IMSCRIBE_P5')),  p5 = getpref(mfilename(),'IMSCRIBE_P5');     end
if (ispref(mfilename(),'IMSCRIBE_P6')),  p6 = getpref(mfilename(),'IMSCRIBE_P6');     end

% --- Check Environment variables ---
if (isempty(p1)), p1 = getenv('IMSCRIBE_P1'); end % location for TEMPLATE Localizer images 
if (isempty(p2)), p2 = getenv('IMSCRIBE_P2'); end % location for TEMPLATE FOV/ROI 
if (isempty(p3)), p3 = getenv('IMSCRIBE_P3'); end % location for TARGET Localizer images 
if (isempty(p4)), p4 = getenv('IMSCRIBE_P4'); end % scratch space folder for temporary files 
if (isempty(p5)), p5 = getenv('IMSCRIBE_P5'); end % RT export folder (for auto load of latest Dicoms) 
if (isempty(p6)), p6 = getenv('IMSCRIBE_P6'); end % Dicom export folder (for generating Dicoms with ROI mask) 

% --- Look for default named folders ---
homepath = getenv('HOME');
is_Windows = ~isempty(findstr(computer(),'PCWIN'));
if (is_Windows), dicompath = 'X:';            diskCpath = 'C:'; 
else,            dicompath = '/mnt/rtexport'; diskCpath = '/mnt/disk_c'; end
if (isempty(p1)), if (isfolder([homepath filesep() 'Data/Imscribe/Templates'])), p1 = [homepath filesep() 'Data/Imscribe/Templates']; end; end
if (isempty(p2)), if (isfolder([homepath filesep() 'Data/Imscribe/Templates'])), p2 = [homepath filesep() 'Data/Imscribe/Templates']; end; end
if (isempty(p3)), if (isfolder([homepath filesep() 'Data/Imscribe/Targets'])),   p3 = [homepath filesep() 'Data/Imscribe/Targets'];   end; end
if (isempty(p4)), if (isfolder([homepath filesep() 'Data/Imscribe/scratch'])),   p4 = [homepath filesep() 'Data/Imscribe/scratch'];   end; end
if (isempty(p5)), if (isfolder([dicompath filesep() 'RTexport_Current'])),       p5 = [dicompath filesep() 'RTexport_Current'];       end; end
if (isempty(p6)), if (isfolder([diskCpath filesep() 'MedCom/temp/CDR_OFFLINE'])), p6 = [diskCpath filesep() 'MedCom/temp/CDR_OFFLINE']; end; end

if (~isdir(p1)), p1 = pwd(); end
if (~isdir(p2)), p2 = pwd(); end
if (~isdir(p3)), p3 = pwd(); end
if (~isdir(p4)), p4 = pwd(); end
%if (~isdir(p5)), p5 = p1; end % these aren't critical and hang if remotely mounted folder is not mounted
%if (~isdir(p6)), p6 = p4; end


% % --- DO NOT EDIT BELOW HERE!!---
% if (~isdir(p4))
%     fprintf(2,'ERROR: The "scratch" folder (%s) does not exist!\n',p4);
%     fprintf(2,'         Edit %s.m to specify default folder locations.\n',mfilename('fullpath'));
%     fprintf(2,'         Alternatively, you can set the following environment variables:\n');p1 = getenv('IMSCRIBE_P1'); % location for TEMPLATE Localizer images 
%     fprintf(2,'           IMSCRIBE_P1  location for TEMPLATE Localizer images\n');
%     fprintf(2,'           IMSCRIBE_P2  location for TEMPLATE FOV/ROI\n');
%     fprintf(2,'           IMSCRIBE_P3  location for TARGET Localizer images\n');
%     fprintf(2,'           IMSCRIBE_P4  scratch space folder for temporary files\n'); 
%     fprintf(2,'           IMSCRIBE_P5  RT export folder (for auto load of latest Dicoms)\n'); 
%     fprintf(2,'           IMSCRIBE_P6  Dicom export folder (for generating Dicoms with ROI mask)\n'); 
%     fprintf(2,'         Restart the Imscribe GUI after fixing this.\n');
%     p4 = '';
%     return
% end

fprintf(1,'Imscribe path1 = "%s"\n',p1);
fprintf(1,'Imscribe path2 = "%s"\n',p2);
fprintf(1,'Imscribe path3 = "%s"\n',p3);
fprintf(1,'Imscribe path4 = "%s"\n',p4);
fprintf(1,'Imscribe path5 = "%s"\n',p5);
fprintf(1,'Imscribe path6 = "%s"\n',p6);

return