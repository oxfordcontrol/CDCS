function X = projectK(X,K,useDual)

%project a vector onto a cone K defined in SEDUMI format

if(nargin < 3)
    useDual = 0;
end

if(useDual)
    dual = 1; primal = 0;
else
    primal = 1; dual = 0;
end

assert(xor(primal,dual));

blockIdx = 1;

% FREE CONE PROJECTION (R^n)
if(isfield(K,'f') && K.f > 0)
    if(primal)
        X{blockIdx} = projectRn(X{blockIdx});
    end
    if(dual)
        X{blockIdx} = projectZero(X{blockIdx});
    end
    blockIdx = blockIdx + 1;
end

% ZERO CONE PROJECTION ({0}^n)
if(isfield(K,'z') && K.z > 0)
    if(primal)
        X{blockIdx} = projectZero(X{blockIdx});
    end
    if(dual)
        X{blockIdx} = projectRn(X{blockIdx});
    end
    blockIdx = blockIdx + 1;
end

% POSITIVE ORTHANT PROJECTION (R^n_+)
if(isfield(K,'l') && K.l > 0)
    X{blockIdx} = projectPosOrthant(X{blockIdx}); %self dual cone
    blockIdx = blockIdx + 1;
end

% SOC PROJECTION
if (isfield(K,'q') && any(K.q))
    for i = (1:length(K.q))
        X{blockIdx} = projectSOC(X{blockIdx});  %self dual cone
        blockIdx = blockIdx + 1;
    end
end

% PSD CONE PROJECTION (S_+^n)
%project everything onto the PSD cone (self-dual)
nPosEigs = 0;
if (isfield(K,'s') && any(K.s))
    for i = (1:length(K.s))
        [X{blockIdx},temp] = projectPSD(X{blockIdx});  %self dual cone
        nPosEigs = nPosEigs + temp;
        blockIdx = blockIdx + 1;
    end
end

% fprintf('Number of positive eigenvalues = %i\n',nPosEigs);

end


%------------------------------------
% ACTUAL PROJECTION FUNCTIONS
%------------------------------------

% ZERO CONE
function X = projectZero(X)
    X(:) = 0;
end

% FREE CONE
function X = projectRn(X)
    return;
end

% POSITIVE ORTHANT
function X = projectPosOrthant(X)
    X = max(X,0);
end

% SECOND ORDER CONE
function X = projectSOC(X)
    s = X(1);
    x = X(2:end);
    normx = norm(x,2);

    if(normx <= s)
        X = X;
    elseif(normx <= -s)
        X = X.*0;
    else
        t = (normx + s)/2;
        X = [t;(t/normx).*x];
    end
end

% POSITIVE SEMIDEFINITE CONE
function [S,nPosEigs] = projectPSD(S)
    %scalar case
    if(numel(S) == 1)
        nPosEigs = double(S > 0);
        S = max(0,S);
        return;
    end
    % matrix case
    S = S./2;
    S = S+S.';
    [U,E] = eig(S);
    ind = diag(E)>0;
    UsE = U(:,ind)*sqrt(E(ind,ind));
    S   = UsE*UsE.';
    nPosEigs = nnz(ind);   
end
