function opts = setUserOpts(opts,userOpts)

% SETUSEROPTS
%
% Set user options

if ~isstruct(userOpts)
    error('Input ''options'' must be a structure.');
else
    fnames = fieldnames(userOpts);
    for n=1:length(fnames)
        if isfield(opts,fnames{n})
            opts.(fnames{n}) = userOpts.(fnames{n});
        else
            warning('Option ''%s'' is unknown and will be ignored.',fnames{n})
        end
    end
end

% make solver lowercase
opts.solver = lower(opts.solver);