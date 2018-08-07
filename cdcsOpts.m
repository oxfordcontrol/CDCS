function options = cdcsOpts

% CDCSOPTS
%
% Default options for <a href="matlab:help('cdcs')">CDCS</a>.
%
% Generic solver options
% ----------------------
% options.solver     = 'hsde';    % which solver (primal/dual/hsde/sos)
% options.relTol     = 1e-4;      % tolerance
% options.rescale    = true;      % scale data to improve convergence
% options.verbose    = 1;         % print or silent
% options.dispIter   = 50;        % print every dispIter iterations
% options.maxIter    = 1000;      % max # iterations
% 
% Chordal decomposition options
% -----------------------------
% options.chordalize = 1;         % how to decompose the constraints (0/1/2)
% options.yPenalty   = true;      % add penalty term for Y block to cost
% options.completion = true;      % complete the unused entries of the decomposed
%                                   primal variable
%
% ADMM penalty options
% --------------------
% options.rho        = 1;         % penalty parameter
% options.adaptive   = true;      % adaptive penalty factor?
% options.tau        = 2;         % increase factor for adaptive penalty scheme 
%                                   (must be > 1)
% options.mu         = 10;        % ratio of residuals for adaptive penalty scheme
% options.rhoMax     = 1e6;       % maximum penalty parameter
% options.rhoMin     = 1e-6;      % minimum penalty parameter
% options.rhoIt      = 10;        % if pres/dres>mu (<mu) mu for rhoIt iterations,
%                                   increase (decrease) options.rho
% 
% Advanced options
% ----------------
% options.KKTfact    = 'blk';     % Options for KKT systems
%                                 %  a) 'blk': block elimination,
%                                 %  b) 'ldl': ldl factor,
%                                 %  c) 'inv': invert
%
% See also CDCS


% Create options structure
options.solver     = 'hsde';    % which solver (primal/dual/hsde)
options.relTol     = 1e-4;      % tolerance
options.rescale    = true;      % scale data to improve convergence
options.verbose    = 1;         % print or silent
options.dispIter   = 50;        % print every dispIter iterations
options.maxIter    = 1000;      % max # iterations
options.chordalize = 1;         % how to decompose the constraints (1/2)
options.yPenalty   = true;      % add penalty term for Y block to cost
options.completion = true;      % complete the unused entries of the primal variable
options.rho        = 1;         % penalty parameter
options.adaptive   = true;      % adaptive penalty factor?
options.mu        = 2;          % increase/decrease factor for adaptive penalty scheme (must be > 1)
options.nu         = 10;        % ratio of residuals for adaptive penalty scheme
options.rhoMax     = 1e6;       % maximum penalty parameter
options.rhoMin     = 1e-6;      % minimum penalty parameter
options.rhoIt      = 10;        % if pres/dres>mu (<mu) mu for rhoIt iterations, adapt rho
options.KKTfact    = 'blk';     % Options for KKT systems

% the following is for self-dual embedding
options.alpha      = 1.8;       % over relaxation, must lie in (0,2); 1.8 is used in SCS

