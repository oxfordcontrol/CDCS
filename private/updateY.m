function [Y,others] = updateY(X,Y,Z,rho,others,projK)

% UPDATEY(X,Y,Z,rho,others,projK)
% Update block Y for sparse conic ADMM solver: a projection on cones

Y.blk = projK(X.blk,Z.blk,rho);
Y.vec = flatten(Y.vec,Y.blk);

end