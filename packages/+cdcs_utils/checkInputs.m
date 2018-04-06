function [At,b,c,K,opts] = checkInputs(At,b,c,K,opts)

%sparsify and vectorize everything
At = sparse(At); 
b = sparse(b(:));

if(isstruct(c))
    c.cx = sparse(c.cx(:));
    c.cz = sparse(c.cz(:));
else
    c = sparse(c(:));
end


%check that only free, zero, non-negative, quadratic cone and SDP variables
%are included
if(~all(ismember(fieldnames(K),{'f','l','q','s'})))
    error('Unsupported cone constraint types.');
end

%basic problem dimensions
[n,m] = size(At);
opts.n_init = n;
opts.m_init = m;

% Set cone
nConeVars = 0;
if(isfield(K,'f') && ~isempty(K.f) && K.f > 0)
    nConeVars = nConeVars + K.f;
else
    K.f = 0;
end

if(isfield(K,'l') && ~isempty(K.l) && K.l > 0)
    nConeVars = nConeVars + K.l;
else
    K.l = 0;
end

if (isfield(K,'q') && ~isempty(K.q) && max(K.q) > 0)
    K.q = K.q(K.q~=0);
    nConeVars = nConeVars + sum(K.q);
else
    K.q = 0;
end

if (isfield(K,'s') && ~isempty(K.s) && max(K.s) > 0)
    K.s = K.s(K.s~=0);
    nConeVars = nConeVars + sum(K.s.^2);
else
    K.s = 0;
end

% Check if cone dimensions match with data
assert(nConeVars == opts.n_init);
