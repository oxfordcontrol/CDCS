function [At,b,c,K,Ech,stuff,opts] = preprocess(At,b,c,K,opts)

% CDCS/packages/+cdcs_sos/PREPROCESS
% rescale data for homogeneous self-dual embedding

% Import functions
import cdcs_utils.svecData
import cdcs_utils.chordalDecomposition

% Any variables at all?
nonSDP = K.f + K.l + sum(K.q);
if nonSDP + sum(K.s) == 0
    error('No variables in your problem?')
end

% no chordal decompositon
Ech = [];
stuff = [];

% svec form
[At,c,~,~] = svecData(At,c,K);
[opts.n,opts.m] = size(At);

% original norms, used repeatly in convergence checking
opts.nAt_init = norm(At,'fro');
opts.nb_init  = norm(b,'fro');
opts.nc_init  = norm(c,'fro');


%--------------------------------------------
% Reorder the PSD cones
%--------------------------------------------
% reorder the variables, such that A = [A1 A2], A2*A2' is diagonal
reorder = 1;
if reorder == 1
    [At,c,K,opts] = ConsReorder(At,c,K,opts);
end

%--------------------------------------------
% Rescale data
%--------------------------------------------
% Rescale
[At,b,c,K,opts] = rescaleData(At,b,c,K,opts); % doesn't work well!


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

% norms, used repeatly in convergence checking
opts.nAt = norm(At,'fro');
opts.nb  = norm(b,'fro');
opts.nc  = norm(c,'fro');



