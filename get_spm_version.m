function [SPMver,spmpath,tissuepath] = get_spm_version(desiredversion)
% Get the version of SPM in your PATH
%
% Syntax:
%   [SPMver,spmpath,tissuepath] = get_spm_version([desiredversion])
% Inputs:
%   desiredversion - specific version of SPM needed (.e.g. 'SPM8')
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

if (nargin < 1), desiredversion = ''; end

if (isdeployed)  % can't check PATH in mcc state, assume it was compiled w/ SPM in path!
    if (~isempty(desiredversion))
        SPMver = desiredversion;
    else
        SPMver = 'SPM8';
    end
    return
end

spmfile = which('spm');
spmpath = fileparts(spmfile);
tissuepath = [spmpath filesep() 'tpm' filesep()];

if isempty(spmfile)
    SPMver = '';
    fprintf(2,'ERROR: Did not find SPM in your MATLAB path\n');
    return
elseif contains(upper(spmpath),'SPM8')
    SPMver = 'SPM8';
elseif contains(upper(spmpath),'SPM12')
    SPMver = 'SPM12';
else
    fprintf(2,'ERROR: Unrecognized version of SPM in "%s"\n',spmfile);
    SPMver = '';
    return
end
if (~isempty(desiredversion))
    if (~contains(desiredversion,SPMver))
        fprintf(2,'ERROR: SPM version in your PATH is %s. It needs to be %s\n',SPMver,desiredversion);
        SPMver = '';
        return
    end
end
end
