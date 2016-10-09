function [Y,others] = updateY(X,Y,Z,rho,others,projK)

% UPDATEY(X,Y,Z,rho,others,projK)
% Update block Y for sparse conic ADMM solver: a projection on cones

% Import functions
import cdcs_utils.flatten

% Operate
Y.blk = projK(X.blk,Z.blk,rho);
Y.vec = cdcs_utils.flatten(Y.vec,Y.blk);

end