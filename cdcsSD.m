function [x,y,info] = cdcsSD(At,b,c,K,userOpts)
% cdcsSD(At,b,c,K,Opts)
% Solve sparse conic program (SDPs in particular) in self-dual embedding.

import cdcs_utils.*
%--------------------------------------------
% Some default parameters
%--------------------------------------------
opts.maxIter    = 1000;      % max # iterations
opts.rescale    = 1;         % scale data or not 
opts.verbose    = 1;         % print or silent
opts.dispIter   = 100;       % print every dispIter iterations
opts.convIter   = 1;         % check for convergence (residual) every this num iterations
opts.relTol     = 1e-3;      % tolerance
opts.rho        = 1;         % penalty parameter
opts.adaptive   = true;      % adaptive penalty factor?  adpative scheme requires the residuals in Byod's paper
opts.tau        = 2;         % increase factor for adaptive penalty scheme (must be > 1)
opts.mu         = 10;        % ratio of residuals for adaptive penalty scheme
opts.rhoMax     = 1e6;       % maximum penalty parameter
opts.rhoMin     = 1e-6;      % minimum penalty parameter
opts.rhoIt      = 5;         % if pres/dres>mu (<mu) mu for rhoIt iterations, adapt rho
opts.alpha      = 1.5;       % over relaxation, must lie in (0,2); 1.5 is recommanded
opts.initVars   = [];        % initial variables
opts.KKTfact    = 'blk';     %'blk': block elimination, 'ldl': ldl factor, 'inv': invert
opts.chordalize = 1;         % how to decompose the constraints (1/2)
opts.completion = true;      % complete the unused entries of the decomposed primal variable
opts.solver     = 'hsde';    % homogeneous self-dual embedding    

if(nargin >= 5)
    opts = setstructfields(opts,userOpts);
end

%% --------------------------------------------
% Setup
%--------------------------------------------
% Print nice welcoming header
myline1 = [repmat('=',1,86),'\n'];
myline2 = [repmat('-',1,86),'\n'];
if opts.verbose
    fprintf(myline1)
    fprintf('CDCS - Homogeneous Self-Dual Embedding v1.0, by Giovanni and Yang \n')
    fprintf('Setting up...')
end

% start timing
proctime = tic;

%% sparsify everything, scale the data, and make chordalization
[At,b,c,K,opts] = checkInputs(At,b,c,K,opts);
[At,b,c,K,opts] = rescaleDataH(At,b,c,K,opts);    % considered the matching matrix H
Kold = K;
[At,b,c,K,Ech,cd,chstuff] = chordalize(At,b,c,K,opts);

% norms of problem data, used in the calculation of convergence repeatedly
opts.nAt = norm(At,'fro');opts.nb  = norm(b,'fro'); opts.nc  = norm(c,'fro');

%% scale the decomposed problem, 
% if using this funtion, please use factorMatrix_scale.m since we have modified matrix H
% [At,b,c,opts] = rescaleDecData(At,b,c,Kold,Ech,opts);  % this function doesn't work well

%% Make operators for ADMM and Initialization
[step1,step2,step3,checkConv] = makeADMM(At,b,c,K,cd,Ech,opts);

[opts.n,opts.m] = size(At);  opts.nXk = length(Ech);
[u,v,X]         = makeVariables(K,opts.initVars,opts);

% Time and display
proctime = toc(proctime);
if opts.verbose
    fprintf('done in %.4f seconds.\n',proctime);
    fprintf(myline1)
    fprintf(' iter |   pres   |   dres   |   pcost   |   dcost   |   gap    |   rho    | time (s) |\n')
    fprintf(myline2)
end

%% --------------------------------------------
%  Iterations of ADMM
%  --------------------------------------------
admmtime = tic;
for iter = 1:opts.maxIter
    
    %% ADMM iterations
    uold = u;
    hatu = step1(u,v,opts.rho);           % Projection to affine set  
    u    = step2(hatu,v,X,opts.rho);      % Projection to cone set
    v    = step3(v,hatu,u,opts.rho);      % update variable v, played as multipliers

    %% check for converengce and detect infeasible problems
    if iter==1 || ~mod(iter,opts.convIter)
        [Flag,pcost,dcost,presi,dresi,gap,opts] = checkConv(u,v,iter,admmtime,opts,hatu,uold); 
        if(Flag.stop)
            break;
        end
    end
end
admmtime = toc(admmtime);

% need a psd completion and scale back to original variables

%% Output information
info.iter      = iter;
info.setupTime = proctime;
info.admmTime  = admmtime;
info.runTime   = proctime+admmtime;
info.pcost     = pcost;
info.dcost     = dcost;
info.pres      = presi;
info.dres      = dresi;
info.gap      = gap;

if Flag.isConverged
    info.feasCode = 'feasible';
    x = u.x/u.tau;
    y = (u.y/u.tau./opts.scaleFactors.E1)./opts.scaleFactors.sc_c;
end

if Flag.pinfeas
   info.feasCode = 'primal infeasible';
   x = NaN;               %% no feasible primal variable exists
   y = u.y/(b'*u.y);      
end

if Flag.dinfeas
   info.feasCode = 'dual infeasible';
   x = u.x/(-c'*u.x);    
   y = NaN;               %% no feasible dual variable exists
end

% GF: max iteration reached?
if iter==opts.maxIter && ~Flag.isConverged
    info.feasCode = 'Max number of iterations reached';
    x = u.x/u.tau;
    y = (u.y/u.tau./opts.scaleFactors.E1)./opts.scaleFactors.sc_c;
end

% Print summary
fprintf(myline1)
fprintf(' SOLUTION SUMMARY:\n')
fprintf('------------------\n')
fprintf(' Number of iterations : %15.d\n',iter)
fprintf([' Problem status       : ',info.feasCode,'\n']);
fprintf(' Primal cost          : %15.4e\n',info.pcost)
fprintf(' Dual cost            : %15.4e\n',info.dcost)
fprintf(' Primal residual      : %15.4e\n',info.pres)
fprintf(' Dual residual        : %15.4e\n',info.dres)
fprintf(' Setup time (s)       : %15.4e\n',proctime)
fprintf(' Admm time  (s)       : %15.4e\n',admmtime)
fprintf(' Total time (s)       : %15.4e\n',info.runTime)
fprintf(myline1)

end
