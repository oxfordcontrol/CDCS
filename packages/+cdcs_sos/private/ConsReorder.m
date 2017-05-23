function [At,c,K,opts] = ConsReorder(At,c,K,opts)
% Reorder the constraints such that A = [A1 A2]. A2*A2' is diagonal

    nCone     = length(K.s);
    nConeVars = cumsum([0+K.f+K.f+K.q, K.s.*(K.s+1)/2]);
    OrthFlag  = zeros(nCone,1);
    OrthNum   = 0;
    for i = 1:nCone  %% focusing on PSD cones
        tmpAt = At(nConeVars(i)+1:nConeVars(i+1),:);
        if isdiag(tmpAt'*tmpAt)   %% nondiagonal part
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

    
    Index = 1;
    ReOrderData = zeros(length(c),1);
    for i = 1:nCone
        Num = nConeVars(ReOrder(i)+1) - nConeVars(ReOrder(i));
        ReOrderData(Index:Index+Num-1)   = nConeVars(ReOrder(i))+1:nConeVars(ReOrder(i)+1);
        Index = Index+Num;
    end
    
    c = c(ReOrderData);
    At = At(ReOrderData,:);
    
    if ~issparse(At) || ~issparse(c)
        c = sparse(c);
        At = sparse(At);
    end
end

