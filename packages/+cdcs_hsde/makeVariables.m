function [u,v,X] = makeVariables(K,initVars,opts)
% MAKEVARIABLES(K,initVars)
% Make initial variables

% GF:
% Initialize with all zeros since we know it is a feasible point.
% Initializing kappa=tau=1 gives annoying NaN in first iteration...

if isempty(initVars)
    %initialize vectorized versions of all variables
    
    u.x  = zeros(opts.n,1)+1;              % free variables 
    u.xh = zeros(opts.nXk,1)+1;            % cone variables: vectorized version
    X    = cdcs_utils.makeConeVariables(K);% cone variables: matrix version
    u.y  = zeros(opts.m,1);                % free variables
    u.v  = zeros(opts.nXk,1);              % free variables
    u.tau = 0;

    v.x  = zeros(opts.n,1);                % similar to u
    v.xh = zeros(opts.nXk,1);              % 
    v.y  = zeros(opts.m,1);                % 
    v.v  = zeros(opts.nXk,1);              % 
    v.kappa = 0;
    
else
    % initialized vectorized variables with used input
    % TO DO LATER
    
end




end