function [X,others] = updateX(X,Y,Z,rho,others,Ech,K,projAffine,opts)

% UPDATEX(X,Y,Z,rho,others,Ech,K,projAffine,opts)
% Update block X for sparse conic ADMM solver: a projection on affine
% constraints

% Import functions
import cdcs_utils.blockify

% Project onto affine constraints
[x,y] = projAffine(Y.vec,Z.vec,rho);

% Blockify 

% Set variables
X.vec(:) = x;
others.dual = y;

%blockify E*x
Ex = x(Ech);%X.vec(Ech);
X.blk = blockify(X.blk,Ex,K);


end