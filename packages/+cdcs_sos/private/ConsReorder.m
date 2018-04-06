function [At,c,K,opts] = ConsReorder(At,c,K,opts)
% Reorder the constraints such that A = [A1 A2]. A2*A2' is diagonal
% Only handle the PSD cone, and leave the free and second-order cone
% unchanged

%% find the diagonal part
    nCone     = length(K.s);
    % GF on 24 July 2017: how about linear variables?
    nConeVars = cumsum([K.f+K.l+K.q, K.s.*(K.s+1)/2]);
    OrthFlag  = zeros(nCone,1);
    OrthNum   = 0;
    for i = 1:nCone  %% focusing on PSD cones
        tmpAt = At(nConeVars(i)+1:nConeVars(i+1),:);
        if isdiag(tmpAt'*tmpAt)   %% diagonal part
            OrthFlag(i) = 1;
            OrthNum = OrthNum+ nConeVars(i+1) - nConeVars(i);
        end
    end
    opts.sos.OrthFlag = OrthFlag;
    opts.sos.Orth = OrthNum;         %% number of orthogonal constraints
    opts.sos.NonOrth = nConeVars(end) - OrthNum; 
    
    [~,ReOrder] = sort(OrthFlag);    %% reorder the variables x
    K.s         = K.s(ReOrder);
    opts.sos.ReOrder = ReOrder;

    %% reorder the data
    Index = 1;
    ReOrderData = zeros(nConeVars(end)-nConeVars(1),1);
    for i = 1:nCone
        Num = nConeVars(ReOrder(i)+1) - nConeVars(ReOrder(i));
        ReOrderData(Index:Index+Num-1)   = nConeVars(ReOrder(i))+1:nConeVars(ReOrder(i)+1);
        Index = Index+Num;
    end
    ReOrderData = [(1:nConeVars(1))';ReOrderData];
    c = c(ReOrderData);
    At = At(ReOrderData,:);
    
    %% sparsify
    if ~issparse(At) || ~issparse(c)
        c = sparse(c);
        At = sparse(At);
    end
end

