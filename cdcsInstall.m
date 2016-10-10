function cdcsInstall

% CDCSINSTALL
%
% Install CDCS and compile required binaries.
%
% See also CDCS, CDCSTEST


% First, compile some files from CSparse
here = pwd;
cd('include')
cs_install
cd(here);
movefile(['include',filesep,'cs_lsolve.mex*'], ...
         [here,filesep,'packages',filesep,'+cdcs_utils']);
movefile(['include',filesep,'cs_ltsolve.mex*'], ...
         [here,filesep,'packages',filesep,'+cdcs_utils']);

% Then compile some mex files from this package
cd(['packages',filesep,'+cdcs_utils',filesep,'private'])
if (~isempty (strfind (computer, '64')))
    mexcmd = 'mex -largeArrayDims' ;
else
    mexcmd = 'mex' ;
end
eval([mexcmd, ' svec.c']);
eval([mexcmd, ' smat.c']);
cd(here)

% Finally add to path and save
addpath([here,filesep,'packages']);
addpath(here);
savepath

fprintf('\nCompilation completed successfully.\n');

end
