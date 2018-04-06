function [At,b,c,K,opts] = rescaleData(At,b,c,K,opts)

% CDCS/packages/+cdcs_sos/RESCALEDATA.m
% Try to rescale data to get nicer convergence properties.
% Assumes no zero columns in At, otherwise will fail.

% Parameters
[n,m] = size(At);
min_scale = 1e-3;
max_scale = 1e+3;
minScaleRowAt = min_scale .*sqrt(m);
maxScaleRowAt = max_scale .*sqrt(m);
minScaleColAt = min_scale .*sqrt(n);
maxScaleColAt = max_scale .*sqrt(n);

if opts.rescale
    
    % Similar scaling strategy to SCS, but scale columns of A (rows of At) to
    % have unit norm. So, flip order of operations compared to SCS!
    
    % Scale rows of At
    % must take mean over cone to preserve cone membership
    % But not true for free and linear blocks: can scale individual variables!
    D = full(sqrt(sum(At.*At,2)));            % norm of cols of A (col vec)
    count = 0;
    if K.f>0
        nvars = K.f;
        % No need for mean in free cone!
        % D(count+1:count+nvars) = mean(D(count+1:count+nvars));
        count = count + nvars;
    end
    if K.l>0
        nvars = K.l;
        % No need for mean in linear cone!
        % D(count+1:count+nvars) = mean(D(count+1:count+nvars));
        count = count + nvars;
    end
    if sum(K.q)>0
        for i = 1:length(K.q)
            nvars = K.q(i);
            D(count+1:count+nvars) = mean(D(count+1:count+nvars));
            count = count + nvars;
        end
    end
    if sum(K.s)>0
        for i = 1:length(K.s)
            % use svec'ed variables
            nvars = K.s(i)*(K.s(i)+1)/2;
            D(count+1:count+nvars) = mean(D(count+1:count+nvars));
            count = count + nvars;
        end
    end
    
    D(D>maxScaleRowAt) = maxScaleRowAt;     % set upper bound
    D(D<minScaleRowAt) = minScaleRowAt;     % set lower bound
    
    if issparse(At)
        At = spdiags(1./D,0,count,count)*At;
    else
        % Old code: with bsxfun
        At = bsxfun(@rdivide,At,D);             % divide row i of At by D(i)
    end
    
    % Scale cols of At
    E = full(sqrt(sum(At.*At,1)));            % norm of rows of A (row vec)
    E(E>maxScaleColAt) = maxScaleColAt;       % set upper bound
    E(E<minScaleColAt) = minScaleColAt;       % set lower bound
    if issparse(At)
        nE = length(E);
        At = At*spdiags(1./E.',0,nE,nE);
    else
        % Old code: with bsxfun
        At = bsxfun(@rdivide,At,E);               % divide col i of At by E(i)
    end
    
    % Find mean row and col norms for scaled A = At.'
    M = At.*At;
    meanRowNormA = mean(sqrt(sum(M,1)));
    meanColNormA = mean(sqrt(sum(M,2)));
    
    % 2) Scale b
    b = b./E(:);
    nm_b = sqrt(b.'*b);
    sc_b = meanColNormA / max(nm_b, min_scale);
    b = b.*sc_b;
    
    % 3) Scale c
    c = c./D;
    nm_c = sqrt(c.'*c);
    sc_c = meanRowNormA / max(nm_c, min_scale);
    c = c.*sc_c;
    
    % 4) Save rescaling variables
    opts.scaleFactors.D = D(:);
    opts.scaleFactors.E = E(:);
    opts.scaleFactors.sc_b = sc_b;
    opts.scaleFactors.sc_c = sc_c;
    opts.scaleFactors.sc_cost = sc_c*sc_b;
    
    
else
    
    % Dummy scale factor
    opts.scaleFactors.D = 1;
    opts.scaleFactors.E = 1;
    opts.scaleFactors.sc_b = 1;
    opts.scaleFactors.sc_c = 1;
    opts.scaleFactors.sc_cost = 1;
    
end

