function [x,y,z,info,opts] = setOutputs(X,Y,Z,K,info,opts)

% SETOUTPUTS.M
%
% Set outputs in sedumi format using positive semidefinite completion algorithms

% Import functions
import cdcs_utils.makeConeVariables
import cdcs_utils.flatten
import cdcs_utils.blockify

x = zeros(opts.n_init,1);
y = zeros(opts.m_init,1);
z = zeros(opts.n_init,1);


%% vector value
xtmp = Y.x/Y.tau;
ytmp = Y.y/Y.tau;
ztmp = Z.x/X.tau; 

% block value
[xmat,zmat] = makeConeVariables(K);        % blockified
xmat  = blockify(xmat,xtmp,K);
zmat  = blockify(zmat,ztmp,K);
xtmp = flatten(xmat,xmat,0); 
ztmp = flatten(zmat,zmat,0); 

% Scale solution
xtmp = (xtmp./opts.scaleFactors.D)./opts.scaleFactors.sc_b;
ytmp = (ytmp./opts.scaleFactors.E)./opts.scaleFactors.sc_c;
ztmp = (ztmp.*opts.scaleFactors.D)./opts.scaleFactors.sc_c;

% Assign solution to used variables
x(opts.usedvars) = xtmp;
z(opts.usedvars) = ztmp;
y = ytmp;

% END FUNCTION
end