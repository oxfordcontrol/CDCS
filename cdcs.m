function [x,y,z,info] = cdcs(At,b,c,K,userOpts,initVars)

% CDCS
%
% Syntax:
%
% [x,y,z,info] = CDCS(At,b,c,K,opts)
%
% Solve a sparse conic program using chordal decomposition for the semidefinite
% cones and ADMM. CDCS solves the primal (P) or dual (D) standard forms of
% the conic problem,
%
%         min <c,x>                             min -<b,y>
%   (P)   s.t. Ax = b,                  (D)     s.t. A^Ty - c = z
%         x \in K                               z \in K*
%
% where A,b and c are the problem date and K is the cone (K* is the dual cone).
% The standard form to be solved is specified by the "solver" field of the
% options structure:
%
%   opts.solver = 'primal' (default): solve the problem in primal standard form
%   opts.solver = 'dual'            : solve the problem in dual standard form
%
% The chordal decomposition can be carried out in two ways, specified by the
% "chordalize" option:
%
%   opts.chordalize = 1 (default): split the data equally between the cliques
%   opts.chordalize = 2          : assign data to one clique only
%
% <a href="matlab:help('cdcsOpts')">Click here for a complete list of options</a>.
%
% The output structure 'info' contains the following information:
%
% info.problem: - 0: CDCS terminated succesfully 
%               - 1: the maximum number of iterations was reached
%               - 2: the ADMM iterations terminated succesfully, but the positive 
%                 matrix completion algorithm threw an error
% info.iter: number of iterations
% info.cost: terminal cost
% info.pres: terminal primal ADMM residual
% info.dres: terminal dual ADMM residual
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
% Solver options
%============================================
tstart = tic;
opts = cdcsOpts;


%============================================
% Setup
%============================================
% Set user options
if(nargin >= 5)
    opts = setUserOpts(opts,userOpts);
end

% Checks on specified solver type and method
if ~any(strcmpi(opts.solver,{'primal','dual'}))
    error('Unknown opts.solver. Please use "primal" or "dual".')
end

% Print nice welcoming header
if opts.verbose
    myline1 = [repmat('=',1,64),'\n'];
    fprintf(myline1)
    fprintf('CDCS by G. Fantuzzi, Y. Zheng -- v1.0\n')
    fprintf(myline1)
    %fprintf('Using method %i for the %s standard form.\n',opts.method,opts.solver)
    fprintf('Initializing CDCS...')
end

% start timing
proctime = tic;

% sparsify everything, check cone constraints, rescale
[At,b,c,K,opts] = checkInputs(At,b,c,K,opts);
[At,b,c,K,opts] = rescaleData(At,b,c,K,opts);
[At,b,c,K,opts] = splitBlocks(At,b,c,K,opts);
[opts.n,opts.m] = size(At);

% chordal decomposition
Kold = K;
[At,b,c,K,Ech,cd,chstuff] = chordalize(At,b,c,K,opts);

% basic decomposed problem dimensions:  no. of cones, no. of vectorized conic
% variables, and no. of free primal variables
opts.p = (K.f>0) + (K.l>0) + length(K.q)*(sum(K.q)>0) + length(K.s)*(sum(K.s)>0);
opts.nXk = length(Ech);
opts.nX  = size(At,1);

% Initial variables
if (nargin < 6); initVars=[]; end
[X,Y,Z,others] = makeVariables(K,initVars,opts);

% Make operators for ADMM
[step1,step2,step3,checkConv] = makeADMM(At,b,c,K,cd,Ech,opts);

% Time setup and display
proctime = toc(proctime);
if opts.verbose
    myline2 = [repmat('-',1,64),'\n'];
    fprintf('done in %.4f seconds.      \n',proctime);
    fprintf('Standard form          : %s\n',opts.solver);
    fprintf('Chordalization method  : %i\n',opts.chordalize);
    fprintf('Adaptive penalty       : %i\n',opts.adaptive);
    fprintf('Scale data             : %i\n',opts.rescale);
    fprintf('Free variables         : %i                \n',K.f);
    fprintf('Non-negative variables : %i                \n',K.l);
    fprintf('Second-order cones     : %i (max. size: %i)\n',length(K.q),max(K.q));
    fprintf('Semidefinite cones     : %i (max. size: %i)\n',length(K.s),max(K.s));
    fprintf('Affine constraints     : %i                \n',opts.m);
    fprintf('Consensus constraints  : %i                \n',sum(accumarray(Ech,1)));
    fprintf(myline1)
end

%============================================
% Run ADMM
%============================================
% Display
if opts.verbose
    fprintf(' iter |   pres   |   dres   |    cost    |   rho    | time (s) |\n')
    fprintf(myline2)
end

admmtime = tic;
opts.feasCode = 1;
for iter = 1:opts.maxIter
    
    % Save current iterate for convergence test
    YOld = Y;
    
    % Update block variables
    [X,others] = step1(X,Y,Z,opts.rho,others);
    [Y,others] = step2(X,Y,Z,opts.rho,others);
    [Z,others] = step3(X,Y,Z,opts.rho,others);
    
    % log errors / check for convergence
    [isConverged,log,opts] = checkConv(X,Y,Z,YOld,others,iter,admmtime,opts);
    if isConverged
        opts.feasCode = 0;
        break;
    end
end
admmtime = toc(admmtime);
if opts.verbose
    fprintf(myline1)
end

%============================================
% Outputs
%============================================
% Variables in sedumi format
posttime = tic;
[x,y,z,opts] = setOutputs(X,Y,Z,others,Kold,cd,Ech,chstuff,opts);
posttime = toc(posttime);

% Info
info.problem = opts.feasCode;              % diagnostic code
info.iter    = iter;                       % # of iterations
info.cost    = log.cost;                   % terminal cost
info.pres    = log.pres;                   % terminal primal ADMM res
info.dres    = log.dres;                   % terminal dual ADMM res
info.time.setup   = proctime;              % setup time
info.time.admm    = admmtime;              % ADMM time
info.time.cleanup = posttime;              % post-processing time
info.time.total   = toc(tstart);           % total CPU time

% Print summary
if opts.verbose
    fprintf(' SOLUTION SUMMARY:\n')
    fprintf('------------------\n')
    fprintf(' Termination code     : %11.1d\n',info.problem)
    fprintf(' Number of iterations : %11.d\n',iter)
    fprintf(' Cost                 : %11.4e\n',info.cost)
    fprintf(' Primal residual      : %11.4e\n',info.pres)
    fprintf(' Dual residual        : %11.4e\n',info.dres)
    fprintf(' Setup time   (s)     : %11.4e\n',proctime)
    fprintf(' ADMM  time   (s)     : %11.4e\n',admmtime)
    fprintf(' Cleanup time (s)     : %11.4e\n',posttime)
    fprintf(' Total time   (s)     : %11.4e\n',info.time.total)
    fprintf(myline1)
end


end











