function [At,b,c,K,opts] = rescaleData(At,b,c,K,Ech,opts)

% +cdcs_hsde/RESCALEDATA.m
% Try to rescale data to get nicer convergence properties.
% taking the matching matrix H into account
% Assumes no zero columns in At, otherwise will fail.

% Parameters
[n,m] = size(At);
min_scale = 1e-3;
max_scale = 1e+3;
% minScaleRowAt = min_scale .*sqrt(m);
% maxScaleRowAt = max_scale .*sqrt(m);
% minScaleColAt = min_scale .*sqrt(n);
% maxScaleColAt = max_scale .*sqrt(n);

minScaleRowAt = min_scale .*sqrt(m+length(Ech));
maxScaleRowAt = max_scale .*sqrt(m+length(Ech));
minScaleColAt = min_scale .*sqrt(n+length(Ech));
maxScaleColAt = max_scale .*sqrt(n+length(Ech));

% -----------------------------------------------------------------
% Matrix before scaling: [At H'] [b] [c]
%                        [0  -I] [0] [0]
%
% Matrix after scaling:  M =[D1AtE1 D1H'E2]  [rho*E1b]  [sigma*D1c]
%                           [0      -D2I2E]  [0      ]  [0        ]
%
% The columns of M and c should have Euclidean norm close to one.
% The rows of M and b should have similar norms.
% -----------------------------------------------------------------

if opts.rescale
    
    % -----------------------------------------------------
    % Scale rows of [At H']
    %               [0  -I]  norms are similar to each other
    % -----------------------------------------------------
    
    % norm of rows of [At H']
    D1 = full(sqrt(sum(At.*At,2) + accumarray(Ech,1)));
    
%    count = 0;
    % No need to take mean 
%     if K.f>0
%         count = count + K.f;
%     end
% 
%     if K.l>0
%         count = count + K.l;
%     end
%     
%     % mean over quaratic cones
%     if sum(K.q)>0
%         for i = 1:length(K.q)
%             nvars = K.q(i);
%             D1(count+1:count+nvars) = mean(D1(count+1:count+nvars));
%             count = count + nvars;
%         end
%     end
%     
%     % mean over SDP cones (with svec!)
%     if sum(K.s)>0
%         for i = 1:length(K.s)
%             nvars = K.s(i)*(K.s(i)+1)/2;
%             D1(count+1:count+nvars) = mean(D1(count+1:count+nvars));
%             count = count + nvars;
%         end
%    end
     
    D1(D1>maxScaleRowAt) = maxScaleRowAt;     % set upper bound
    D1(D1<minScaleRowAt) = 1;                 % set lower bound
    if issparse(At)
        nD = length(D1);
        At = spdiags(1./D1,0,nD,nD)*At;
    else
        % Old code: with bsxfun
        At = bsxfun(@rdivide,At,D1);              % divide row i of At by D(i)
    end
    % norms of rows of [0 -I] 
    D2 = ones(length(Ech),1);       %% always maintain the membership of cones
    
    % -----------------------------------------------------
    % Scale cols of [D1At D1H']  [c]
    %               [0    -D2I]  [0]  norms close to one
    % -----------------------------------------------------
    
    % norm of cols of [At] 
    E1 = full(sqrt(sum(At.*At,1)));            
    E1(E1>maxScaleColAt) = maxScaleColAt;       % set upper bound
    E1(E1<minScaleColAt) = 1;                   % set lower bound
    if issparse(At)
        nE = length(E1);
        At = At*spdiags(1./E1.',0,nE,nE);
    else
        At = bsxfun(@rdivide,At,E1);                % divide col i of At by E(i)
    end
    
    % cols of [D1H']
    %         [-D2I]
    % rows of [HD1 -ID2]
    tmp = (1./D1).^2;
    E2 = sqrt(tmp(Ech)+1); 
    %E2 = ones(length(Ech));
    
    % --------------------------------------------------------
    % Find mean row and col norms for scaled  [D1AtE1 D1H'E2] = [D1AtE1 D1H'E2]
    %                                         [0      -D2IE2]   [0       -IE2 ]
    % --------------------------------------------------------
    M1 = At.*At;
    M2 = sum(M1,2) + accumarray( Ech,(1./E2).^2 )./(D1.^2);
    tmp = (1./D1).^2;
    meanRowNormAt = mean([sqrt(M2); 1./E2]);
    meanColNormAt = mean([sqrt(sum(M1,1)),sqrt(tmp(Ech)./(E2.^2)+(1./E2).^2)']);
    
    % Scale b
    b = b./E1(:);
    nm_b = sqrt(b.'*b);
    sc_b = meanRowNormAt / max(nm_b, min_scale);
    b = b.*sc_b;
    
    % Scale c
    c = c./D1;
    nm_c = sqrt(c.'*c);
    sc_c = meanColNormAt / max(nm_c, min_scale);
    c = c.*sc_c;
    
    % --------------------------------------------------------
    % Save rescaling variables
    % --------------------------------------------------------
    % D2 and E2 do not appear in the computation later? 
    % Try to incoporate the influence of D1, D2, E1, E2
    opts.scaleFactors.D1 = 1./D1(:);
    opts.scaleFactors.E1 = 1./E1(:);
    opts.scaleFactors.D2 = 1./D2(:);  %% identity vector
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