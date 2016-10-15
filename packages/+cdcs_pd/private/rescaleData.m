function [At,b,c,K,opts] = rescaleData(At,b,c,K,opts)

% CDCS/packages/+cdcs_pd/RESCALEDATA.m
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
    
    % Similar scaling strategy to SCS
    
    % Scale rows of At 
    % must take mean over cone to preserve cone membership
    D = full(sqrt(sum(At.*At,2)));            % norm of cols of A (col vec)
    count = 0;
    if K.f>0
        nvars = K.f;
        D(count+1:count+nvars) = mean(D(count+1:count+nvars));
        count = count + nvars;
    end
    if K.l>0
        nvars = K.l;
        D(count+1:count+nvars) = mean(D(count+1:count+nvars));
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
            nvars = K.s(i)^2;
            D(count+1:count+nvars) = mean(D(count+1:count+nvars));
            count = count + nvars;
        end
    end
    
    D(D>maxScaleRowAt) = maxScaleRowAt;     % set upper bound
    D(D<minScaleRowAt) = 1;                 % set lower bound
    At = bsxfun(@rdivide,At,D);             % divide row i of At by D(i)
    
    % Scale cols of At
    E = full(sqrt(sum(At.*At,1)));            % norm of rows of A (row vec)
    E(E>maxScaleColAt) = maxScaleColAt;       % set upper bound
    E(E<minScaleColAt) = 1;                   % set lower bound
    At = bsxfun(@rdivide,At,E);               % divide col i of At by E(i)
    
    
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

