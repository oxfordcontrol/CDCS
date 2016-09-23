function admmPDCPInstall

%ADMMPDCPINSTALL
%
% Install admmPDCP and compile required binaries.


% First, compile some files from CSparse
here = pwd;
cd('include')
cs_install
cd(here);
movefile(['include',filesep,'cs_lsolve.mex*'],[here,filesep,'private']);
movefile(['include',filesep,'cs_ltsolve.mex*'],[here,filesep,'private']);

% Then compile some mex files from this package
cd('private')
if (~isempty (strfind (computer, '64')))
    mexcmd = 'mex -largeArrayDims' ;
else
    mexcmd = 'mex' ;
end
eval([mexcmd, ' svec.c']);
eval([mexcmd, ' smat.c']);
cd(here)

% Finally add to path and save
addpath(here);
savepath

fprintf('\nCompilation completed successfully.\n');

end
