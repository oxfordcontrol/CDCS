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


% Decompose cost vector
if opts.chordalize == 0
    % Do nothing!
    cd = c;
elseif opts.chordalize == 1
    % Decompose equally
    IA  = accumarray(E,1);
    cd  = c./IA; cd  = cd(E);
elseif opts.chordalize == 2
    % Decompose using only last entry
    nv  = length(E);
    cd  = zeros(nv,1);
    [U,IA] = unique(E,'last');
    cd(IA,:) = c(U,:);
else
    error('Unknown chordal decomposition method.')
end

% Create some useful variables. xsvec is initialized so it is obvious which
% entries have not been computed if psdCompletion is not called upon user
% request. This is also good for compatibility with YALMIP: used variables are
% set to NaN!
xtemp = zeros(opts.n,1);
ztemp = zeros(opts.n,1);
xsvec = NaN(chstuff.totvars,1);            % svec format
zsvec = zeros(chstuff.totvars,1);          % svec format
[xmat,zmat] = makeConeVariables(K);        % blockified

% Extract variables, depending on algorithm used
if strcmpi(opts.solver,'primal')
    xsvec(chstuff.usedvars) = X.vec;
    y(:) = others.dual;
    zsvec(chstuff.usedvars) = accumarray(E,Z.vec+cd.*opts.yPenalty);
    zsvec = zsvec./(1+opts.yPenalty);
    zmat = blockify(zmat,zsvec,K);
    ztemp = flatten(ztemp,zmat,0);
    
elseif strcmpi(opts.solver,'dual')
    xsvec(chstuff.usedvars) = others.dual;
    y(:) = X.vec(1:opts.m);
    zsvec(chstuff.usedvars) = accumarray(E,Y.vec);
    zmat = blockify(zmat,zsvec,K);
    ztemp = flatten(ztemp,zmat,0);
    
end

% Positive semidefinite completion of x variable
% Only complete if problem successfully solved!
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

% Scale solution
xtemp = (xtemp./opts.scaleFactors.D)./opts.scaleFactors.sc_b;
y = (y./opts.scaleFactors.E)./opts.scaleFactors.sc_c;
ztemp = (ztemp.*opts.scaleFactors.D)./opts.scaleFactors.sc_c;

% Assign solution to used variables
x(opts.usedvars) = xtemp;
z(opts.usedvars) = ztemp;


% END FUNCTION
end