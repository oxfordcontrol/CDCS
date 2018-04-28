function [cliques,Ech,Jch,usedvars,s] = ConeDecomposition(K,opts)


% CDCS/packages/+cdcs_hocp/ConeDecomposition.m
%
% Compute chordal decomposition of SDP cones.

import cdcs_utils.ConeSplit

%--------------------------------------------
% Non SDP variables
%--------------------------------------------
nonSDP =  K.f + K.l + sum(K.q);
nonSDPind = (1:nonSDP)';


% %--------------------------------------------
% % Sparsity pattern, choral extension and maximal cliques
% %--------------------------------------------
% SP = spones(spones(Cs) + sparse(sum(spones(Ats),2)));  % vector of 1s and 0s

% Decompose each block
nCones  = length(K.s);
count   = 0;
countE  = 0;
usedvars = [];
s       = cell(nCones,1);
cliques = cell(nCones,1);
Ech     = {};
Jch     = {};
for i = 1:nCones
    
    % make block
    %         n = K.s(i)^2;
    %         Blk = reshape(SP(count+1:count+n),K.s(i),K.s(i));
    n = K.s(i)*(K.s(i)+1)/2;
%    Blk = smat(SP(count+1:count+n));
    
    % Cone decomposition
    if isfield(opts,'partition') && ~isempty(opts.partition)
        alpha = opts.partition;
%     elseif isfield(opts,'NoP') && ~isempty(opts.NoP)
%         alpha = ones(9,1)*floor(K.s(i)/opts.NoP);
%         alpha = [alpha;K.s(i) - sum(alpha)];
    else
        if K.s(i) > 10
            alpha = ones(opts.NoP-1,1)*floor(K.s(i)/opts.NoP);
            alpha = [alpha;K.s(i) - sum(alpha)];
        else
            alpha = K.s(i);
        end
    end
    cliques{i} = ConeSplit(alpha);
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
        E{j}   = LOCB + nonSDP + countE;
        ind{j} = ind{j} + nonSDP + count;
    end
    
    Ech = [Ech; E];
    Jch = [Jch; ind];
    count = count + n;
    countE = countE + length(nnzind);
end



%--------------------------------------------
% Concatenate mapping indices
%--------------------------------------------
Ech = [nonSDPind; vertcat(Ech{:})];
Jch = [nonSDPind; vertcat(Jch{:})];
s   = horzcat(s{:});
usedvars = [nonSDPind; usedvars];