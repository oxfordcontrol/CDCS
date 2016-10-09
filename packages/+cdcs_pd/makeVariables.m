function [X,Y,Z,others] = makeVariables(K,initVars,opts)

% MAKEVARIABLES(initVars,K)
% Make initial variables
%
% X.vec: vectorized variable for affine constraints
% Y.vec: list of vectorized sub-variables for conic projections
% Z.vec: vectorized Lagrange multipliers
% X.blk: blockified variable for affine constraints - must match Y.blk
% Y.blk: blockified subvariables, ready for projection onto cones
% Z.blk: blockified Lagrange multipliers


if isempty(initVars)
    %initialize vectorized versions of all variables
    if strcmpi(opts.solver,'primal')
        X.vec = sparse(opts.nX,1);
        others.dual = sparse(opts.m,1);
    elseif strcmpi(opts.solver,'dual')
        X.vec = sparse(opts.m+opts.nXk,1);
        others.dual = sparse(opts.nX,1);
    end
    Y.vec = sparse(opts.nXk,1);
    Z.vec = sparse(opts.nXk,1);
    
    %initialize blockified versions of all variables
    [X.blk,Y.blk,Z.blk] = cdcs_utils.makeConeVariables(K);
    
else
    % initialized vectorized variables with used input
    X.vec = initVars.X.vec;
    Y.vec = initVars.Y.vec;
    Z.vec = initVars.Z.vec;
    X.blk = initVars.X.blk;
    Y.blk = initVars.Y.blk;
    Z.blk = initVars.Z.blk;
    
end




end