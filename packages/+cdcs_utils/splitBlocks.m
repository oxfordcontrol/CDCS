function [At,b,c,K,opts] = splitBlocks(At,b,c,K,opts)

% Try to split semidefinite blocks into smaller connected components
% Return updated data with split cones and a vector of indices opts.sort such
% that
%
% x_new = x_old(opts.sort)
%
% Given x_new, the original variable x_old can be reconstructed with
%
%  >> x_old = zeros(<appropriate size>,1);
%  >> x_old(opts.sort) = x_new;

% Some variables
cnt = K.f + K.l + sum(K.q);
nBlk = length(K.s);

% Sparsity pattern of semidefinite blocks
csdp = c(cnt+1:end);
Atsdp = At(cnt+1:end,:);
SP = spones(spones(csdp) + sparse(sum(spones(Atsdp),2)));

% Loop over blocks
R = []; rowcnt = 0;
C = []; colcnt = 0;
s = [];
for i = 1:nBlk
    
    % Find connected components of block
    m = K.s(i);
    B = reshape(SP(rowcnt+1:rowcnt+m^2),m,m);
    [tags,nComponents] = connectedComponents(B);
    
    % If more than one component, split.
    if nComponents>1
               
        % Find component dimensions, permutation matrix and indices for block
        % pattern
        p = zeros(m,1); 
        count = 0; 
        I = []; J = [];
        for j = 1:nComponents
            % Dimension and permutation vector
            ind = find(tags==j);
            dim = numel(ind);
            s = [s, dim];
            lind = count+1:count+dim;
            p(lind) = ind;
            
            % Find indices for block pattern
            [row,col] = meshgrid(lind);
            I = [I; row(:)]; 
            J = [J; col(:)];
            
            % Update counter
            count = count+dim;
            
        end
        
        
        % Project onto connected components
        P = sparse(p,1:m,1,m,m);                    % permutation of cone
        B = sparse(I,J,1,m,m);                      % permuted connected components
        e = find(B); ne = length(e);
        E = sparse(e,1:ne,1,m^2,ne);  % projection onto component blocks
        [Ri,Ci] = find(kron(P,P)*E);
        
        % Store indices
        R = [R; Ri+rowcnt];
        C = [C; Ci+colcnt];      
        
        % Update counters
        rowcnt = rowcnt + m^2;
        colcnt = colcnt + ne;
        
    else
        % Nothing to do
        s = [s, m];
        Ri = (1:m^2)';
        R = [R; Ri+rowcnt];
        C = [C; Ri+colcnt];
        rowcnt = rowcnt + m^2;
        colcnt = colcnt + m^2;
    end
    
end


% Update problem data
M = sparse(C,R,1,colcnt,rowcnt);    % since work with transpose At!
At = [At(1:cnt,:); M*Atsdp];
c  = [c(1:cnt); M*csdp];
K.s = s;

% Find used variables
opts.usedvars = [(1:cnt)'; R+cnt];



% END MAIN
end








% ============================================================================ %
%                               NESTED FUNCTION                                %
% ============================================================================ %

function [tags,nComponents] = connectedComponents(B)

    % Find connected components of sparsity pattern matrix B

    L = size(B,1); % number of vertex

    % Breadth-first search:
    tags = zeros(1,L); % all vertex unexplored at the begining
    rts = [];
    ccc = 0; % connected components counter
    while true
        ind = find(tags==0);
        if ~isempty(ind)
            fue  = ind(1); % first unexplored vertex
            rts  = [rts fue];
            list = [fue];
            ccc  = ccc+1;
            tags(fue) = ccc;
            while true
                list_new = [];
                for lc=1:length(list)
                    p   = list(lc); % point
                    cp  = find(B(p,:)); % points connected to p
                    cp1 = cp(tags(cp)==0); % get only unexplored vertecies
                    tags(cp1) = ccc;
                    list_new  =[list_new cp1];
                end
                list = list_new;
                if isempty(list)
                    break;
                end
            end
        else
            break;
        end
    end
    
    % set number of components
    nComponents = max(tags);
end