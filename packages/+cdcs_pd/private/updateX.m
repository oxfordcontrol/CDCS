function [X,others] = updateX(X,Y,Z,rho,others,Ech,K,projAffine,opts)

% UPDATEX(X,Y,Z,rho,others,Ech,K,projAffine,opts)
% Update block X for sparse conic ADMM solver: a projection on affine
% constraints

% Import functions
import cdcs_utils.blockify

% Project onto affine constraints
[x,y] = projAffine(Y.vec,Z.vec,rho);

% Blockify - primal and dual options
if strcmpi(opts.solver,'primal')
    
    % Set variables
    X.vec(:) = x;
    others.dual = y;
    
    %blockify E*x
    Ex = x(Ech);%X.vec(Ech);
    X.blk = blockify(X.blk,Ex,K);
    
elseif strcmpi(opts.solver,'dual')
    
    % Set variables
    vk = Y.vec + (Z.vec - x(Ech))./rho; 
    X.vec(:) = [y; vk];
    others.dual = x;
    
    % Blockify vk inside X.blk to use in updateY
    X.blk = blockify(X.blk,vk,K);
    
end

end