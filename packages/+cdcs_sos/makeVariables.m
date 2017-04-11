function [X,Y,Z,others] = makeVariables(K,initVars,opts)
% [u,v,X] = makeVariables(K,initVars,opts)
% uhat -- X; u --> Y, v --> Z;
% MAKEVARIABLES(K,initVars)
% Make initial variables


if isempty(initVars)
    %initialize vectorized versions of all variables
    X.x   = zeros(opts.n,1);                % cone variables 
    X.y   = zeros(opts.m,1);                % free variables
    X.tau = sqrt(opts.n+opts.m+1);          % set as sqrt(n+m+1) in SCS?
    
    
    Y.x   = zeros(opts.n,1);                % cone variables 
    Y.y   = zeros(opts.m,1);                % free variables
    Y.tau = sqrt(opts.n+opts.m+1);

    Z.x     = zeros(opts.n,1);                % similar to Y
    Z.y     = zeros(opts.m,1);                % 
    Z.kappa = 0;
    
    others.X    = cdcs_utils.makeConeVariables(K);% cone variables: matrix version

else
    % initialized vectorized variables with used input
    % TO DO LATER
    
end

end