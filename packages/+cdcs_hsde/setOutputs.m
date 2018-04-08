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
y(:) = (Y.y./Y.tau.*opts.scaleFactors.E1)./opts.scaleFactors.sc_c;
zsvec(chstuff.usedvars) = accumarray(E,Y.v./Y.tau);
zsvec(chstuff.usedvars) = (zsvec(chstuff.usedvars)./opts.scaleFactors.D1)./opts.scaleFactors.sc_c;
zmat = blockify(zmat,zsvec,K);
ztemp = flatten(ztemp,zmat,0);

% Positive semidefinite completion of x variable
% Only complete if problem successfully solved!

tmpScale = max(opts.scaleFactors.D1)./opts.scaleFactors.sc_b;    % keep the unscaled variables for completion
xsvec(chstuff.usedvars) = (Y.x./Y.tau.*opts.scaleFactors.D1)./opts.scaleFactors.sc_b./tmpScale;
xmat  = blockify(xmat,xsvec,K);
xtemp = flatten(xtemp,xmat,0);                   % in sedumi format for psdCompletion
if info.problem==0 && opts.completion==1 && opts.chordalize~=0
    try
        % This will give an error if one PSD cone can in fact be split into multiple
        % separate cones
        xtemp = psdCompletion(xtemp,K,chstuff.cliques); % psdComplete!
    catch
        warning('CDCS:psdCompletion',...
            ['Aborting matrix completion algorithm due to a problem.\n'...
            'Variables in the positive semidefinite cones will be ',...
            'returned without completion.']);
        info.problem = 4;
    end
elseif info.problem==3
    warning('CDCS:psdCompletion',...
        ['CDCS reached the maximum number of iterations, and will not attempt\n',...
        'to complete the positive semidefinite variable. Your output will most\n',...
        'likely contain NaNs!']);
end

% scale back the solutions


% Assign solution to used variables
x(opts.usedvars) = xtemp.*tmpScale;
z(opts.usedvars) = ztemp;


%% check individual residual for the consensus constraint Hkx = x or zk = vk 
    % added by Yang 09/03/2018
    p = length(others.X);   % number of cliques
    Xktemp = blockify(others.X,Y.x(E),others.K);
    Xhtemp = blockify(others.X,Y.xh,others.K);
    Zktemp = blockify(others.X,Z.xh,others.K);
    Vktemp = blockify(others.X,Y.v,others.K);
    cpres = zeros(p,2);
    cdres = zeros(p,2);
    pEig  = zeros(p,1);
    dEig  = zeros(p,1);
    for k = 1:p
        cpres(k,1) = norm(Xktemp{k} - Xhtemp{k},'fro');
        cpres(k,2) = cpres(k,1)./ max([norm(Xktemp{k},'fro'), norm(Xhtemp{k},'fro'),1]);
        
        cdres(k,1) = norm(Zktemp{k} - Vktemp{k},'fro');
        cdres(k,2) = cdres(k,1)./ max([norm(Zktemp{k},'fro'), norm(Vktemp{k},'fro'),1]);
        
        pEig(k) = min(eig(Xktemp{k}));
        dEig(k) = min(eig(Vktemp{k}));
    end
    info.cpres = max(cpres); % maximal primal consensus residual;
    info.cdres = max(cdres); % maximal dual consensus residual;
    
    info.pEig = min(pEig);
    info.dEig = min(dEig);


% END FUNCTION
end