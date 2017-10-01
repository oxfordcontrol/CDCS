function [At,b,c,K,Ech,stuff,opts] = preprocess(At,b,c,K,opts)

% CDCS/packages/+cdcs_hsde/PREPROCESS
% Chordal decomposition of SDP constraints & rescale data for homogeneous
% self-dual embedding

% Import functions
import cdcs_utils.svecData
import cdcs_utils.chordalDecomposition

% Any variables at all?
nonSDP = K.f + K.l + sum(K.q);
if nonSDP + sum(K.s) == 0
    error('No variables in your problem?')
end

%--------------------------------------------
% Chordal decomposition
%--------------------------------------------

% Returns the submatrix Ats, Cs of data for PSD variables with svec
% instead of vec, which is the standard input. Also return modified
% matrices At,c that account for svec operation.
[At,c,Ats,Cs] = svecData(At,c,K);
totvars = size(At,1);

% If required and have SDP variables, decompose
if ~isempty(K.s) && any(K.s~=0) && opts.chordalize~=0
    [cliques,Ech,Jch,usedvars,s] = chordalDecomposition(Ats,Cs,K);
    
else
    cliques = [];
    Ech = (1:totvars)';
    Jch = (1:totvars)';
    usedvars = (1:totvars)';
    s = K.s;   
end

% original norms, used repeatly in convergence checking
opts.nAt_init = norm(At,'fro');
opts.nb_init  = norm(b,'fro');
opts.nc_init  = norm(c,'fro');

%--------------------------------------------
% Rescale data
%--------------------------------------------
[At,b,c,K,opts] = rescaleData(At,b,c,K,Jch,opts);


%--------------------------------------------
% Set new At, c, K.s (remove unused variables)
%--------------------------------------------
At  = At(usedvars,:);
c   = c(usedvars);
K.s = s;

if opts.rescale == 1
    opts.scaleFactors.D1 = opts.scaleFactors.D1(usedvars);
end

% Check if At,b,C are indeed sparse - if not, make full for speed!
[n,m] = size(At);
densityTol = 0.6;
if nnz(At)/(m*n) > densityTol
    At = full(At);
end
if nnz(b)/m > densityTol
    b = full(b);
end
if nnz(c)/n > densityTol
    c = full(c);
end


%--------------------------------------------
% Set stuff
%--------------------------------------------
stuff.cliques  = cliques;
stuff.usedvars = usedvars;
stuff.totvars  = totvars;

% norms, used repeatly in convergence checking
opts.nAt = norm(At,'fro');
opts.nb  = norm(b,'fro');
opts.nc  = norm(c,'fro');
