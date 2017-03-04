function [x,y,z,info,opts] = setOutputs(X,Y,Z,K,info,opts)

% SETOUTPUTS.M
%
% Set outputs in sedumi format using positive semidefinite completion algorithms

% Import functions
import cdcs_utils.makeConeVariables
import cdcs_utils.flatten
import cdcs_utils.blockify


%% vector value
x = Y.x/Y.tau;
y = Y.y/Y.tau;
z = Z.x/X.tau; 

% block value
[xmat,zmat] = makeConeVariables(K);        % blockified
xmat  = blockify(xmat,x,K);
zmat  = blockify(zmat,z,K);
x = flatten(xmat,xmat,0); 
z = flatten(zmat,zmat,0); 

% END FUNCTION
end