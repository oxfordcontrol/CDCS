function cdcsTest

% CDCSTEST
%
% Run some example to test CDCS.
%
% See also CDCS

% Preliminaries
clc;
opts.maxIter = 1e+3;
opts.relTol  = 1e-3;


% ---------------------------------------------------------------------------- %
%           SDP with block-arrow sparsity pattern (multiple cones)
% ---------------------------------------------------------------------------- %
% Parameters
m = 100;                        % # constraints
nCones = 1;                     % # cones with block-arrow sparsity pattern
nBlk = [15, 20];                % # diagonal blocks for each cone
BlkSize = [10, 5];              % block size for each cone
ArrowHead = [5, 5];             % arrow head size for each PSD cone

% Setup
fprintf('\nSetting up random block-arrow SDP, m=%i...',m);
tsetup = tic;
[At,b,c,K] = blockArrowMultCones(m,nCones,nBlk,BlkSize,ArrowHead);
tsetup = toc(tsetup);
fprintf('done in %.2f seconds. \n',tsetup);

% solution by admm
opts.solver = 'primal';
cdcs(At,b,c,K,opts);
opts.solver = 'dual';
cdcs(At,b,c,K,opts);
% opts.solver = 'hsde'; % note: current implementation of hsde fails; need debug
% cdcs(At,b,c,K,opts);


% ---------------------------------------------------------------------------- %
%                       Conic program with banded SDP
% ---------------------------------------------------------------------------- %
% Parameters
m   = 300;                      % # constraints
K.f = 23;                       % # free variables
K.l = 150;                      % # non-negative variables
K.q = [15, 30];                 % # second-order cones
K.s = [75 33];                    % PSD cones
bandWidth = [10 5];              % bandWidth for SDP cones

% Setup
fprintf('\nSetting up random conic problem with banded SDP cones, m=%i...',m);
tsetup = tic;
[At,b,c,K] = bandedSDP(m,K,bandWidth);
tsetup = toc(tsetup);
fprintf('done in %.2f seconds. \n',tsetup);

% solution by admm
opts.solver = 'primal';
cdcs(At,b,c,K,opts);
opts.solver = 'dual';
cdcs(At,b,c,K,opts);
opts.solver = 'hsde';
cdcs(At,b,c,K,opts);

% ---------------------------------------------------------------------------- %
%                              SDPs in SDPLIB
% ---------------------------------------------------------------------------- %
% qap9
fprintf('Testing the SDPLIB problem qap9 \n');
load(['examples',filesep,'qap9.mat'])
opts.solver = 'primal';
cdcs(At,b,c,K,opts);
opts.solver = 'dual';
cdcs(At,b,c,K,opts);
opts.solver = 'hsde';
cdcs(At,b,c,K,opts);

% mcp250-1
fprintf('\nTesting the SDPLIB problem mcp250-1 \n');
load(['examples',filesep,'mcp250-1.mat'])
opts.solver = 'primal';
cdcs(At,b,c,K,opts);
opts.solver = 'dual';
cdcs(At,b,c,K,opts);
opts.solver = 'hsde';
cdcs(At,b,c,K,opts);

% ---------------------------------------------------------------------------- %
%                                   END
% ---------------------------------------------------------------------------- %
fprintf('\n\nCDCS was successfully tested.\n\n')
end



% ============================================================================ %
%                               NESTED FUNCTIONS                               %
% ============================================================================ %

% -------------------
% blockArrowMultCones
% -------------------
function [At,b,c,K] = blockArrowMultCones(m,nCones,nBlk,BlkSize,ArrowHead)

% Setup problem for block-arrow sdp with multiple cones. Inputs:
% m        : number of equality constraints
% nCones   : number of cones with same block-arrow sparsity pattern
% nBlk     : number of diagonal blocks in block-arrow matrix data (vector of length nCones)
% BlkSize  : size of each diagonal block (vector of length nCones)
% ArrowHead: size of head of arrow pattern (vector of length nCones)

% cone
K.f = 0;
K.l = 0;
K.q = 0;
K.s = zeros(1,nCones);

% Sparsity pattern of each cone
Spa = cell(nCones,1);
for k = 1:nCones
    n = nBlk(k)*BlkSize(k)+ArrowHead(k);
    Spa{k} = zeros(n);
    for i = 1:nBlk(k)
        Spa{k}((i-1)*BlkSize(k)+1:BlkSize(k)*i, ...
               (i-1)*BlkSize(k)+1:BlkSize(k)*i) = ones(BlkSize(k));
    end
    Spa{k}(nBlk(k)*BlkSize(k)+1:n,:) = 1;
    Spa{k}(:,nBlk(k)*BlkSize(k)+1:n) = 1;
    
    % Set cone size
    K.s(k) = n;
end

% Data
At = [];
for i = 1:m
    Ai = [];
    for k = 1:nCones
        M = 100*sprandsym(Spa{k});  % random symmetric data with given sparsity pattern
        Ai = [Ai; M(:)];        % concatenate
    end
    At = [At, Ai(:)];
end


% stictly feasible primal point
X = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
    X{k} = Temp(:);
end
b = At'*vertcat(X{:});

% stictly feasible dual point
y = rand(m,1);
S = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
    S{k} = Temp(:);
end
c = vertcat(S{:}) + At*y;

end



% -------------------
% bandedSDP
% -------------------
function [At,b,c,K] = bandedSDP(m,K,bandWidth)

% Setup problem for block-arrow sdp with multiple cones. Inputs:
% m        : number of equality constraints
% K        : the cone
% bandWidth: the bandWidth for each SDP cone (vector of the same length as K.s)

import cdcs_utils.makeConeVariables
import cdcs_utils.projectK
import cdcs_utils.clean
density = 0.35;

% Check cone
if(~all(ismember(fieldnames(K),{'f','l','q','s'})))
    error('Unsupported cone constraint types.');
end
if ~isfield(K,'f') || isempty(K.f)
    K.f = 0;
end
if ~isfield(K,'l') || isempty(K.l)
    K.l = 0;
end
if ~isfield(K,'q') || isempty(K.q)
    K.q = 0;
end
if ~isfield(K,'s') || isempty(K.s)
    K.s = 0;
end

% Check bandwidth
if any(bandWidth > K.s - 1)
    error('Specified bandwidth must be smaller than the size of the SDP cone.')
end

% Data
At = [];
for i = 1:m
    Ai = sprand(K.f,1,density);                 % free vars
    Ai = [Ai; sprand(K.l,1,density)];           % non-negative orthant
    Ai = [Ai; sprand(sum(K.q),1,density)];      % second-order cone
    for k = 1:length(K.s)                       % PSD cone
        n = 2*bandWidth(k) + 1;
        B = 10*rand(K.s(k),n)-50;
        M = spdiags(B,-bandWidth(k):bandWidth(k),K.s(k),K.s(k));    % random data with 
        M = M + M.';                                                % given bandwidth  
        Ai = [Ai; M(:)];                                      % concatenate
    end
    
    % Concatenate
    At = [At, Ai(:)];
end

% stictly feasible primal point
X = makeConeVariables(K);
for i = 1:length(X)
    X{i} = 100.*rand(size(X{i},1),size(X{i},2))-50;
end
X = projectK(X,K,0);         % point in the primal cone
X = cellfun(@(x)x(:),X,'UniformOutput',false);
x = vertcat(X{:});
b = At'*x;

% strictly feasible dual point
y = rand(m,1);
S = makeConeVariables(K);
shift = 0;
for i = 1:length(S)-length(K.s)
    S{i} = 100.*rand(size(S{i},1),size(S{i},2))-50;
    shift = shift + 1;
end
for i = 1:length(K.s)
    n = 2*bandWidth(i) + 1;
    B = rand(K.s(i),n)-1;
    M = spdiags(B,-bandWidth(i):bandWidth(i),K.s(i),K.s(i));    % random data with 
    M = full(M + M.');                                          % given bandwidth 
    S{shift+i} = M + (-min(eig(M))+1)*eye(K.s(i));
end
S = projectK(S,K,1);         % point in the dual cone (with given sparsity pattern)
S = cellfun(@(x)x(:),S,'UniformOutput',false);
s = vertcat(S{:});
c = s + At*y;
c = clean(c,1e-10);

end



% ============================================================================ %
%                        END OF NESTED FUNCTIONS                               %
% ============================================================================ %
