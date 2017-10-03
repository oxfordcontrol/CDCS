function [X,Y,Z,others] = makeVariables(K,initVars,opts)
% [u,v,X] = makeVariables(K,initVars,opts)
% uhat -- X; u --> Y, v --> Z;
% MAKEVARIABLES(K,initVars)
% Make initial variables

% GF:
% Initialize with all zeros since we know it is a feasible point.
% Initializing kappa=tau=1 gives annoying NaN in first iteration...


if isempty(initVars)
    %initialize vectorized versions of all variables
    X.x  = zeros(opts.nX,1)+1;             % free variables 
    X.xh = zeros(opts.nXk,1)+1;            % cone variables: vectorized version
    X.y  = zeros(opts.m,1);                % free variables
    X.v  = zeros(opts.nXk,1);              % free variables
    X.tau = 0;
    
    Y.x  = zeros(opts.nX,1)+1;             % free variables 
    Y.xh = zeros(opts.nXk,1)+1;            % cone variables: vectorized version
    Y.y  = zeros(opts.m,1);                % free variables
    Y.v  = zeros(opts.nXk,1);              % free variables
    Y.tau = 0;

    Z.x  = zeros(opts.nX,1);               % similar to Y
    Z.xh = zeros(opts.nXk,1);              % 
    Z.y  = zeros(opts.m,1);                % 
    Z.v  = zeros(opts.nXk,1);              % 
    Z.kappa = 0;
    
    others.X    = cdcs_utils.makeConeVariables(K);% cone variables: matrix version

else
    % initialized vectorized variables with used input
    % TO DO LATER
    
end

end