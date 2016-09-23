function opts = admmPDCPopts

% List of default options for admmPDCP
%
% Generic solver options
% ----------------------
% opts.solver     = 'primal';  % which solver (primal/dual)?
% opts.relTol     = 1e-4;      % tolerance
% opts.rescale    = true;      % scale data to improve convergence
% opts.verbose    = 1;         % print or silent
% opts.dispIter   = 50;        % print every dispIter iterations
% opts.maxIter    = 1000;      % max # iterations
% 
% Chordal decomposition options
% -----------------------------
% opts.chordalize = 1;         % how to decompose the constraints (1/2)
% opts.yPenalty   = true;      % add penalty term for Y block to cost
% opts.completion = true;      % complete the unused entries of the decomposed
%                                primal variable
%
% ADMM penalty options
% --------------------
% opts.rho        = 1;         % penalty parameter
% opts.adaptive   = true;      % adaptive penalty factor?
% opts.tau        = 2;         % increase factor for adaptive penalty scheme (must be > 1)
% opts.mu         = 10;        % ratio of residuals for adaptive penalty scheme
% opts.rhoMax     = 1e6;       % maximum penalty parameter
% opts.rhoMin     = 1e-6;      % minimum penalty parameter
% opts.rhoIt      = 10;        % if pres/dres>mu (<mu) mu for rhoIt iterations, adapt rho
% 
% Advanced options
% ----------------
% opts.KKTfact    = 'blk';     % Options for KKT systems
%                              %  a) 'blk': block elimination,
%                              %  b) 'ldl': ldl factor,
%                              %  c) 'inv': invert
%
% See also admmPDCP


% Create options structure
opts.solver     = 'primal';  % which solver (primal/dual)?
opts.relTol     = 1e-4;      % tolerance
opts.rescale    = true;      % scale data to improve convergence
opts.verbose    = 1;         % print or silent
opts.dispIter   = 50;        % print every dispIter iterations
opts.maxIter    = 1000;      % max # iterations
opts.chordalize = 1;         % how to decompose the constraints (1/2)
opts.yPenalty   = true;      % add penalty term for Y block to cost
opts.completion = true;      % complete the unused entries of the primal variable
opts.rho        = 1;         % penalty parameter
opts.adaptive   = true;      % adaptive penalty factor?
opts.tau        = 2;         % increase factor for adaptive penalty scheme (must be > 1)
opts.mu         = 10;        % ratio of residuals for adaptive penalty scheme
opts.rhoMax     = 1e6;       % maximum penalty parameter
opts.rhoMin     = 1e-6;      % minimum penalty parameter
opts.rhoIt      = 10;        % if pres/dres>mu (<mu) mu for rhoIt iterations, adapt rho
opts.KKTfact    = 'blk';     % Options for KKT systems
