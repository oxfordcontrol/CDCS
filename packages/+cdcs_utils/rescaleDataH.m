function [At,b,c,K,opts] = rescaleDataH(At,b,c,K,opts)

% RESCALEDATA.m
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
    
    %% constructing matching matrix H
    Ech = chordalH(At,c,K);
      
    % Scale rows of [At H']
    %               [0  -I]
    % must take mean over cone to preserve cone membership
    D = full(sqrt(sum(At.*At,2)+accumarray(Ech,1)));            % norm of rows of [At H']
    count = 0;
    if K.f>0  % do not need to take mean over free cone
        nvars = K.f;
        count = count + nvars;
    end
    if K.l>0  % do not need to take mean over non-negative cone
        nvars = K.l;
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
    
    D1 = D;%(1:n);
    D2 = ones(length(Ech),1);               % norms of rows of [0 -I]
    
    % Scale cols of [D1At D1H']
    %               [0   -D2I ]
    E = full(sqrt(sum(At.*At,1)));            % norm of cols of [At] (row vec)
    E(E>maxScaleColAt) = maxScaleColAt;       % set upper bound
    E(E<minScaleColAt) = 1;                   % set lower bound
    At = bsxfun(@rdivide,At,E);               % divide col i of At by E(i)
    
    E1 = E; %(1:n);
    E2 = sqrt(ones(length(Ech),1)+(1./D1(Ech)).^2); % rows of [D1H -D2I]
    
    % Find mean row and col norms for scaled  [D1AtE1 D1H'E2]
    %                                         [0      -D2IE2]
    M = At.*At;
    meanRowNormA = 1;%mean(sqrt(sum(M,1))); This number is usually one;
    meanColNormA = mean([sqrt( sum(M,2) + accumarray( Ech,(1./E2).^2 )./(D1.^2));abs(1./E2./D2)]);
    
    % 2) Scale b
    b = b./E1(:);
    nm_b = sqrt(b.'*b);
    sc_b = meanColNormA / max(nm_b, min_scale);
    b = b.*sc_b;
    
    % 3) Scale c
    c = c./D1;
    nm_c = sqrt(c.'*c);
    sc_c = meanRowNormA / max(nm_c, min_scale);
    c = c.*sc_c;
    
    % 4) Save rescaling variables
    opts.scaleFactors.D1 = 1./D1(:);   
    opts.scaleFactors.E1 = 1./E1(:);   
    opts.scaleFactors.D2 = 1./D2(:);   % D2 and E2 do not appear in the computation later !!
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
end

%% nested function

function Ech = chordalH(At,c,K)
% Constructing the matching matrix H

    % Any variables at all?
    if K.f + K.l + sum(K.q) + sum(K.s) == 0
        error('No variables in your problem?')
    end

    % Sort out variables not in SDP cones
    nonSDP = K.f + K.l + sum(K.q);
    nonSDPind = (1:nonSDP)';

    % If have SDP variables to decompose
    if ~isempty(K.s) && any(K.s~=0)

        Cs = c(nonSDP+1:end);Ats = At(nonSDP+1:end,:);         % SDP data
        SP = spones(spones(Cs) + sparse(sum(spones(Ats),2)));  % vector of 1s and 0s

        % Decompose each block
        nCones  = length(K.s);
        countE  = nonSDP;
        countC   = 0;
        cliques = cell(nCones,1);
        Ech     = [];
        for i = 1:nCones

            % make block
            n = K.s(i)^2;
            Blk = reshape(SP(countC+1:countC+n),K.s(i),K.s(i));

            % Chordal decomposition
            cliques{i} = cliquesFromSpMatD(Blk);
            P = cliques{i}.NoC;
            MCO = sparse(K.s(i),P);
            for k = 1:P
                MCO(:,k) = sparse(cliques{i}.Set{k},1,1,K.s(i),1);
            end

            %% coordinates from local copies to global variables, and vice versa
            E = cell(P,1);
            ind = cell(P,1);
            for j = 1:P
                Position = find(MCO(:,j) == 1);
                [cols,rows] = meshgrid(Position,Position);
                E{j} = sub2ind([K.s(i),K.s(i)],rows(:),cols(:)) + countE;   % indices in global variable (including zero entries)
            end


            %%  stack in a vertor form
            E  = full(vertcat(E{:}));   %% \sum n_k^2 * 1, vector
            Ech = [Ech; E];
            countC = countC + n;
            countE = countE + n;
        end

        %--------------------------------------------
        % Concatenate mapping indices
        %--------------------------------------------
        Ech = [nonSDPind; Ech];

    else
        Ech = nonSDPind;    
    end
end



