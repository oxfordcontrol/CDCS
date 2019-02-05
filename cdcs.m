function [x,y,z,info] = cdcs(At,b,c,K,userOpts,initVars)
% CDCS
%
% Syntax:
%
% [x,y,z,info] = CDCS(At,b,c,K,options)
%
% Solve a sparse conic program using chordal decomposition for the positive semidefinite
% cones and ADMM. CDCS solves the primal (P) or dual (D) standard forms of
% the conic problem,
%
%         min <c,x>                             max <b,y>
%   (P)   s.t. Ax = b,                  (D)     s.t. A^Ty + z = c
%         x \in K                               z \in K*
%
% where A, b and c are the problem data and K is the cone (K* is the dual cone).
% CDCS supports the following cones: Free, Linear, second-order,
% Semi-definite, called K.f, K.l, K.q, and K.s.
%
% The standard form to be solved is specified by the "solver" field of the
% options structure:
%
%   options.solver = 'hsde' (default): solve the problem in homogeneous self-dual embedding form
%   options.solver = 'primal'        : solve the problem in primal standard form
%   options.solver = 'dual'          : solve the problem in dual standard form
%   options.solver = 'sos'           : solve the problem arising from Sum-of-squares programs
%
% The chordal decomposition can be carried out in two ways, specified by the
% "chordalize" option:
%
%   options.chordalize = 1 (default): split the data equally between the cliques
%   options.chordalize = 2          : assign data to one clique only
%
% <a href="matlab:help('cdcsOpts')">Click here for a complete list of options</a>.
%
% The output structure 'info' contains the following information:
%
% info.problem: - 0: CDCS terminated successfully 
%               - 1: primal infeasibility detected
%               - 2: dual infeasibility detected
%               - 3: maximum number of iterations reached
%               - 4: the ADMM iterations terminated successfully, but the positive 
%                    matrix completion algorithm threw an error
% info.iter: number of iterations
% info.cost: terminal cost
% info.pres: terminal primal ADMM residual
% info.dres: terminal dual ADMM residual
% info.log : history log of the ADMM residuals, cost, etc.
% info.time: some timing information (setup, ADMM iterations, cleanup, total)
%
% See also CDCSOPTS


% UNDOCUMENTED OPTION:
% An initial guess for the variables used in the decomposed problem can be
% specified using the input "initVars". The variables should be specified in the
% same format as the internal variables - check ./private/makeVariables.m.


% Copyright: G. Fantuzzi [1], Y. Zheng [2], P. Goulart [2],
%            A. Papachristodoulou [2], A. Wynn [1],
%            13 September 2016
%
% [1] Department of Aeronautics, Imperial College London, South Kensington
%     Campus, SW7 2AZ, London, UK.
% [2] Department of Engineering Science, University of Oxford, Parks Road,
%     OX1 3PJ, Oxford, UK


%============================================
% Solver options & import cdcs_utils
%============================================
tstart = tic;
opts = cdcsOpts;
import cdcs_utils.*

%============================================
% Setup
%============================================
% Set user options
if(nargin >= 5)
    opts = setUserOpts(opts,userOpts);
end

% Checks on specified solver type and method
if ~any(strcmpi(opts.solver,{'primal','dual','hsde','sos'}))
    error('Unknown opts.solver. Please use "primal", "dual", "hsde" or "sos".')
end

% Print nice welcoming header
if opts.verbose
    [header,myline1,myline2] = printHeader(opts);
    fprintf(myline1)
    fprintf('CDCS by G. Fantuzzi, Y. Zheng -- v1.0\n')
    fprintf(myline1)
    fprintf('Initializing CDCS...')
end

% start timing
proctime = tic;

% sparsify everything, check cone constraints
[At,b,c,K,opts] = checkInputs(At,b,c,K,opts);
[At,b,c,K,opts] = splitBlocks(At,b,c,K,opts);
[opts.n,opts.m] = size(At);

% rescaling & chordal decomposition for primal/dual/hsde
Kold = K;
[At,b,c,K,Ech,chstuff,opts] = preprocess(At,b,c,K,opts);

% basic decomposed problem dimensions:  no. of cones, no. of vectorized conic
% variables, and no. of free primal variables
opts.p = (K.f>0) + (K.l>0) + length(K.q)*(sum(K.q)>0) + length(K.s)*(sum(K.s)>0);
opts.nXk = length(Ech);
opts.nX  = size(At,1);

% Initial variables
if (nargin < 6); initVars=[]; end
[X,Y,Z,others] = makeVariables(K,initVars,opts);

% Make operators for ADMM
[updateX,updateY,updateZ,checkConvergence] = makeADMM(At,b,c,K,Ech,opts);

% Time setup and display
proctime = toc(proctime);
if opts.verbose
    % Set method to display
    if strcmpi(opts.solver,'hsde')
        method = 'homogeneous self-dual embedding';
    else
        method = opts.solver;
    end
    fprintf('done in %.4f seconds.      \n',proctime);
    fprintf('Algorithm              : %s\n',method);
    %fprintf('Chordalization method  : %i\n',opts.chordalize);
    if any(strcmpi(opts.solver,{'primal','dual'}))
    fprintf('Adaptive penalty       : %i\n',opts.adaptive);
    end
    fprintf('Scale data             : %i\n',opts.rescale);
    fprintf('Free variables         : %i                \n',K.f);
    fprintf('Non-negative variables : %i                \n',K.l);
    fprintf('Second-order cones     : %i (max. size: %i)\n',length(find(K.q ~=0)),max(K.q));
    fprintf('Semidefinite cones     : %i (max. size: %i)\n',length(find(K.s ~=0)),max(K.s));
    fprintf('Affine constraints     : %i                \n',opts.m);
    if any(strcmpi(opts.solver,{'primal','dual','hsde'}))
    fprintf('Consensus constraints  : %i                \n',sum(accumarray(Ech,1)));
    else
    fprintf('Nonorthogonal dimension: %i                \n',opts.sos.NonOrth);
    end
    fprintf(myline1);
    fprintf(header);
    fprintf(myline2);
end

%============================================
% Run ADMM
%============================================
subTime  = zeros(opts.maxIter,3);  % linear proj., conic proj., dual update
log.cost = zeros(opts.maxIter,1);
log.pres = zeros(opts.maxIter,1);
log.dres = zeros(opts.maxIter,1);

admmtime = tic;
for iter = 1:opts.maxIter
    % Save current iterate for convergence test
    YOld = Y;

    % Update block variables
    linearProj = tic;
    [X,others] = updateX(X,Y,Z,opts.rho,others);
    subTime(iter,1) = toc(linearProj);

    conicProj  = tic;
    [Y,others] = updateY(X,Y,Z,opts.rho,others);
    subTime(iter,2) = toc(conicProj);

    dualUpdate  = tic;
    [Z,others]  = updateZ(X,Y,Z,opts.rho,others);
    subTime(iter,3) = toc(dualUpdate);

    % log errors / check for convergence
    [stop,info,log,opts] = checkConvergence(X,Y,Z,YOld,others,iter,admmtime,opts,log);
    if stop
        break;
    end
end
admmtime = toc(admmtime);

%============================================
% Outputs
%============================================
% Variables in sedumi format
posttime = tic;
[x,y,z,info,opts] = setOutputs(X,Y,Z,others,Kold,c,Ech,chstuff,info,opts);
posttime = toc(posttime);

% Info
info.iter    = iter;                       % # of iterations
info.cost    = log.cost(iter);             % terminal cost
info.pres    = log.pres(iter);             % terminal primal ADMM res
info.dres    = log.dres(iter);             % terminal dual ADMM res
info.log.pres     = log.pres(1:iter);      % log of residuals etc
info.log.dres     = log.dres(1:iter);
info.log.cost     = log.cost(1:iter);
if any(strcmpi(opts.solver,{'hsde','sos'}))
    info.log.gap     = log.gap(1:iter);
end
info.time.setup   = proctime;              % setup time
info.time.admm    = admmtime;              % ADMM time
info.time.cleanup = posttime;              % post-processing time
info.time.total   = toc(tstart);           % total CPU time
info.time.subiter = sum(subTime);          % time for each subiteration

% Print summary
if opts.verbose
    fprintf(myline1)
    fprintf(' SOLUTION SUMMARY:\n')
    fprintf('------------------\n')
    fprintf(' Termination code     : %11.1d\n',info.problem)
    fprintf(' Number of iterations : %11.d\n',iter)
    fprintf(' Cost                 : %11.4e\n',info.cost)
    fprintf(' Primal residual      : %11.4e\n',info.pres)
    fprintf(' Dual residual        : %11.4e\n',info.dres)
    fprintf(' Setup time   (s)     : %11.4e\n',proctime)
    fprintf(' ADMM  time   (s)     : %11.4e\n',admmtime)
    fprintf(' Avg. conic proj (s)  : %11.4e\n',info.time.subiter(2)./iter)
    fprintf(' Avg. affine proj (s) : %11.4e\n',info.time.subiter(1)./iter)
    fprintf(' Cleanup time (s)     : %11.4e\n',posttime)
    fprintf(' Total time   (s)     : %11.4e\n',info.time.total)
    fprintf(myline1)
end


end











