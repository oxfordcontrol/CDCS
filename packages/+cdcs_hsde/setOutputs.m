function [x,y,z,info,opts] = setOutputs(X,Y,Z,others,K,c,E,chstuff,info,opts)

% SETOUTPUTS.M
%
% Set outputs in sedumi format using positive semidefinite completion algorithms

% Import functions
import cdcs_utils.makeConeVariables
import cdcs_utils.flatten
import cdcs_utils.blockify
import cdcs_utils.psdCompletion

% Intialize outputs with correct dimensions
x = zeros(opts.n_init,1);
y = zeros(opts.m_init,1);
z = zeros(opts.n_init,1);


% Create some useful variables. xsvec is initialized so it is obvious which
% entries have not been computed if psdCompletion is not called upon user
% request. This is also good for compatibility with YALMIP: used variables are
% set to NaN!
xtemp = zeros(opts.n,1);
ztemp = zeros(opts.n,1);
xsvec = NaN(chstuff.totvars,1);            % svec format
zsvec = zeros(chstuff.totvars,1);          % svec format
[xmat,zmat] = makeConeVariables(K);        % blockified

% Extract variables & scale solution back
xsvec(chstuff.usedvars) = (Y.x./Y.tau.*opts.scaleFactors.D1)./opts.scaleFactors.sc_b;
y(:) = (Y.y./Y.tau.*opts.scaleFactors.E1)./opts.scaleFactors.sc_c;
zsvec(chstuff.usedvars) = accumarray(E,Z.xh./Y.tau);
zsvec(chstuff.usedvars) = (zsvec(chstuff.usedvars).*opts.scaleFactors.D1)./opts.scaleFactors.sc_c;

zmat = blockify(zmat,zsvec,K);
ztemp = flatten(ztemp,zmat,0);
% Positive semidefinite completion of x variable
xmat  = blockify(xmat,xsvec,K);
xtemp = flatten(xtemp,xmat,0);                   % in sedumi format for psdCompletion
if (info.problem==0 || info.problem==3)&& opts.completion==1
    try
        % This will give an error if one PSD cone can in fact be split into multiple
        % separate cones
        xtemp = psdCompletion(xtemp,K,chstuff.cliques); % psdComplete!
    catch
        warning('CDCS:psdCompletion',...
            ['Aborting matrix completion algorithm due to a problem.\n'...
            'Variables in the positive semidefinite cones will be ',...
            'returned without completion.']);
        info.problem = 2;
    end
end

% Assign solution to used variables
x(opts.usedvars) = xtemp;
z(opts.usedvars) = ztemp;


% END FUNCTION
end