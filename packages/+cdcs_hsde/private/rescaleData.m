function [At,b,c,K,opts] = rescaleData(At,b,c,K,Ech,opts)

% +cdcs_hsde/RESCALEDATA.m
% Try to rescale data to get nicer convergence properties.
% taking the matching matrix H into account
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
    
    % -----------------------------------------------------
    % Scale rows of [At H']
    %               [0  -I]
    % -----------------------------------------------------
    
    % norm of rows of [At H']
    D1 = full(sqrt( sum(At.*At,2) + accumarray(Ech,1) ));
    count = 0;
    
    % No need to take mean for K.f and K.l
    if K.f>0
        count = count + K.f;
    end
    
    if K.l>0
        count = count + K.l;
    end
    
    % mean over quaratic cones
    if sum(K.q)>0
        for i = 1:length(K.q)
            nvars = K.q(i);
            D1(count+1:count+nvars) = mean(D1(count+1:count+nvars));
            count = count + nvars;
        end
    end
    
    % mean over SDP cones (with svec!)
    if sum(K.s)>0
        for i = 1:length(K.s)
            nvars = K.s(i)*(K.s(i)+1)/2;
            D1(count+1:count+nvars) = mean(D1(count+1:count+nvars));
            count = count + nvars;
        end
    end
    
    D1(D1>maxScaleRowAt) = maxScaleRowAt;     % set upper bound
    D1(D1<minScaleRowAt) = 1;                 % set lower bound
    At = bsxfun(@rdivide,At,D1);             % divide row i of At by D(i)
    
    % norms of rows of [0 -I]
    D2 = ones(length(Ech),1);
    
    
    % -----------------------------------------------------
    % Scale cols of [D1At D1H']
    %               [0   -D2I ]
    % -----------------------------------------------------
    
    % norm of cols of [At] (row vec)
    E1 = full(sqrt(sum(At.*At,1)));            
    E1(E1>maxScaleColAt) = maxScaleColAt;       % set upper bound
    E1(E1<minScaleColAt) = 1;                   % set lower bound
    At = bsxfun(@rdivide,At,E1);               % divide col i of At by E(i)
    
    % rows of [D1H -D2I]
    E2 = sqrt(ones(length(Ech),1)+(1./D1(Ech)).^2); 
    
    % --------------------------------------------------------
    % Find mean row and col norms for scaled  [D1AtE1 D1H'E2]
    %                                         [0      -D2IE2]
    % --------------------------------------------------------
    M1 = At.*At;
    M2 = sum(M1,2) + accumarray( Ech,(1./E2).^2 )./(D1.^2);
    meanRowNormA = mean(sqrt(sum(M1,1)));
    meanColNormA = mean([sqrt(M2); abs(1./E2./D2)]);
    
    % Scale b
    b = b./E1(:);
    nm_b = sqrt(b.'*b);
    sc_b = meanColNormA / max(nm_b, min_scale);
    b = b.*sc_b;
    
    % Scale c
    c = c./D1;
    nm_c = sqrt(c.'*c);
    sc_c = meanRowNormA / max(nm_c, min_scale);
    c = c.*sc_c;
    
    % --------------------------------------------------------
    % Save rescaling variables
    % --------------------------------------------------------
    % D2 and E2 do not appear in the computation later?
    opts.scaleFactors.D1 = 1./D1(:);
    opts.scaleFactors.E1 = 1./E1(:);
    opts.scaleFactors.D2 = 1./D2(:);
    opts.scaleFactors.E2 = 1./E2(:);
    opts.scaleFactors.sc_b = sc_b;
    opts.scaleFactors.sc_c = sc_c;
    opts.scaleFactors.sc_cost = sc_c*sc_b;
    
    
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

% End function
end