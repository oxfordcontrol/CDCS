function [At,b,c,K,Ech,cd,stuff] = chordalize(At,b,c,K,opts)

% Chordal decomposition of SDP constraints

% Any variables at all?
if K.f + K.l + sum(K.q) + sum(K.s) == 0
    error('No variables in your problem?')
end

% Sort out variables not in SDP cones
nonSDP = K.f + K.l + sum(K.q);
nonSDPind = (1:nonSDP)';


% If have SDP variables to decompose
if ~isempty(K.s) && any(K.s~=0)
    
    %--------------------------------------------
    % Separate part of At, C for PSD vars
    %--------------------------------------------
    % Returns the submatrix Ats, Cs of data for PSD variables with svec
    % instead of vec, which is the standard input. Also return modified
    % matrices At,c that account for svec operation.
    [At,c,Ats,Cs] = svecData(At,c,K);
    totvars = size(At,1);
    
    %--------------------------------------------
    % Sparsity pattern, choral extension and maximal cliques
    %--------------------------------------------
    SP = spones(spones(Cs) + sparse(sum(spones(Ats),2)));  % vector of 1s and 0s
    
    % Decompose each block
    nCones  = length(K.s);
    count   = 0;
    countE  = 0;
    usedvars = [];
    s       = cell(nCones,1);
    cliques = cell(nCones,1);
    Ech     = {};
    for i = 1:nCones
        
        % make block
%         n = K.s(i)^2;
%         Blk = reshape(SP(count+1:count+n),K.s(i),K.s(i));
        n = K.s(i)*(K.s(i)+1)/2;
        Blk = smat(SP(count+1:count+n));
        
        % Chordal decomposition
        cliques{i} = cliquesFromSpMatD(Blk);
        P = cliques{i}.NoC;
        Blk = spones(cliques{i}.idxMatrix+cliques{i}.idxMatrix');   %% chordal sparsity
        ri = cliques{i}.Elem;
        ci = zeros(length(ri),1);
        tmp = 0;
        for k = 1:P
            ci( tmp+1 : tmp+cliques{i}.NoElem(k) )  = k;
            tmp = tmp + cliques{i}.NoElem(k);
        end
        MCO = sparse(ri,ci,1,K.s(i),P);
        
        % Find variables to keep
        nnzind = find(svec(Blk));
        usedvars = [usedvars; nnzind+nonSDP+count];
        
                
        % Indexing to extract local submatrices & split cone
        E = cell(P,1);
        ind = cell(P,1);
        s{i} = zeros(1,P);
        for j = 1:P
            
            %indices in global variable X (including zero entries)
            Position = find(MCO(:,j));
            s{i}(j) = length(Position);             % size of clique cone   
            Position = repmat(Position,1,s{i}(j));
            rw = Position(:);
            cl = Position'; cl = cl(:);
            lti = rw>=cl;                        % lower triangular indices
            rw = rw(lti);
            cl = cl(lti);
            ind{j} = rw+(K.s(i)-cl./2).*(cl-1); % linear indices in svec(X)
            
            % Which entry of list of nonzero elements in global variable?
            [LIA,LOCB] = ismember(ind{j},nnzind);
            E{j} = LOCB + nonSDP + countE;
            
        end
        
        Ech = [Ech; E];
        count = count + n;
        countE = countE + length(nnzind);
    end
   
    
    %--------------------------------------------
    % Remove rows from At, C  corresponding to variables that can be dropped
    % according to the aggregate sparsity pattern in SP
    %--------------------------------------------
    K.s      = horzcat(s{:});
    usedvars = [(1:nonSDP)';usedvars];
    At       = At(usedvars,:);
    c        = c(usedvars);
    
    % Check if At,b,C are indeed sparse - if not, make full for speed!
    [n,m]=size(At);
    densityTol = 0.6;
    if nnz(At)/(m*n) > densityTol
        At = full(At);
    end
    if nnz(b)/m > densityTol
        b = full(b);
    end
    if nnz(c)/n > densityTol
        c = full(c);
    end
    
    %--------------------------------------------
    % Concatenate mapping indices
    %--------------------------------------------
    Ech = [nonSDPind; vertcat(Ech{:})];

    %--------------------------------------------
    % Decompose problem data if necessary
    %--------------------------------------------
    
    if opts.chordalize == 1
        % Decompose equally
        IA  = accumarray(Ech,1);
        cd  = c./IA; cd  = cd(Ech);
        
    elseif opts.chordalize == 2
        % Decompose using only last entry
        nv  = length(Ech);
        cd  = zeros(nv,1);
        [U,IA] = unique(Ech,'last');
        cd(IA,:) = c(U,:);
        
    else
        error('Unknown chordal decomposition method.')
        
    end
    
    
else
    Ech = nonSDPind;
    usedvars = nonSDPind;
    totvars = size(At,1);
    
end

%--------------------------------------------
% Set stuff
%--------------------------------------------
stuff.cliques  = cliques;
stuff.usedvars = usedvars;
stuff.totvars  = totvars;
end