function [At,b,c,K] = bandedSDP(m,K,bandWidth)

% Setup problem for block-arrow sdp with multiple cones. Inputs:
% m        : number of equality constraints
% K        : the cone
% bandWidth: the bandWidth for each SDP cone (vector of the same length as K.s)

import cdcs_utils.makeConeVariables
import cdcs_utils.projectK
import cdcs_utils.clean
density = 1;

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


