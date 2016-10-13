function [At,b,c,opts] = rescaleDecData(At,b,c,K,Ech,opts)

% RESCALEDATA.m
% Try to rescale data to get nicer convergence properties.
% Assumes no zero columns in At, otherwise will fail.
% Considered the matching matrix H, and also need to consider the effects
% of moidifying H when solving the inner system

% Parameters
[n,m] = size(At);
min_scale = 1e-3;
max_scale = 1e+3;
minScaleRowAt = min_scale .*sqrt(m);
maxScaleRowAt = max_scale .*sqrt(m);
minScaleColAt = min_scale .*sqrt(n);
maxScaleColAt = max_scale .*sqrt(n);

A = At';   %% using this format
if opts.rescale
    

    %% construct matrix A0 = [A   ]
    %                        [H -I]
    ns = length(Ech);
    H = sparse([1:ns],Ech,1);
    I = speye(ns);
    A0 = [A, zeros(m,ns);
          H, -I];
    A0 = sparse(A0);  
    %% Scale cols of [A   ]
    %                [H -I], which should have similar norms with b
    % must take mean over cone to preserve cone membership
    
    %D1 = full(sqrt(sum(A.*A,1) + accumarray(Ech,1)'));            % norm of cols of [A;H]
    D = full(sqrt(sum(A0.*A0,1))); 
    D(1:n) = mean(D(1:n));
    D(D>maxScaleRowAt) = maxScaleRowAt;     % set upper bound
    D(D<minScaleRowAt) = 1;                 % set lower bound
    A0 = bsxfun(@rdivide,A0,D);             % divide col i of At by D1(i)
    
    D1 = D(1:n);
    D2 = D(n+1:end);
    %D2 = ones(length(Ech),1);                 % norms of rows of I
    
    %% Scale rows of [A   ]
    %                [H -I]
    E = full(sqrt(sum(A0.*A0,2)));              % norm of row vec of A0
    E(E>maxScaleColAt) = maxScaleColAt;         % set upper bound
    E(E<minScaleColAt) = 1;                     % set lower bound
    A0 = bsxfun(@rdivide,A0,E);                 % divide row i of A0 by E(i)
       
    E1 = E(1:m); 
    E2 = E(m+1:end);
    % H = E2*H*D1 --> THIS SHOULD BE INCLUDED IN SOLVING THE INNER SYSTEM
    
    %% Find mean row and col norms for scaled A = At.'
    M = A0.*A0;
    meanRowNormA = mean(sqrt(sum(M,2))); 
    meanColNormA = mean(sqrt(sum(M,1)));

    %% 2) Scale b
    b = b./E1(:);
    nm_b = sqrt(b.'*b);
    sc_b = meanColNormA / max(nm_b, min_scale);
    b = b.*sc_b;
    
    %% 3) Scale c
    c = c./D1(:);
    nm_c = sqrt(c.'*c);
    sc_c = meanRowNormA / max(nm_c, min_scale);
    c = c.*sc_c;
    
    %% 4) Save rescaling variables
    opts.scaleFactors.D1 = 1./D1(:);   
    opts.scaleFactors.E1 = 1./E1(:);   
    opts.scaleFactors.D2 = 1./D2(:);   
    opts.scaleFactors.E2 = 1./E2(:);   
    opts.scaleFactors.sc_b = sc_b;
    opts.scaleFactors.sc_c = sc_c;
    opts.scaleFactors.sc_cost = sc_c*sc_b;
    
    A = A0(1:m,1:n);
    At = A';
      
else
    
    % Dummy scale factor
    opts.scaleFactors.D1 = 1;
    opts.scaleFactors.E1 = 1;
    opts.scaleFactors.D2 = 1;
    opts.scaleFactors.E2 = 1;
    opts.scaleFactors.sc_b = 1;
    opts.scaleFactors.sc_c = 1;
    opts.scaleFactors.sc_cost = 1;
    
end

