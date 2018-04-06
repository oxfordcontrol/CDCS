function [eta,solInner] = factorMatrix(At,b,c,K,flag)
% generate projector onto affine constraints
% check orthogonality A = [A1 A2], A2 correponds to PSD cone variabels

if strcmpi(flag,'blk')
    [A1,~,D] = diviConstraint(At,K);
    Ddiag = diag(D);   %% store the diagonal elements
    P = 1./(1+Ddiag);

    % First factor the matrix (I+A1'PA1) where P = inv(I+D) is a diagonal
    % matrix; if A1 is empty (all the constraints are orthogonal), then return P
    factors = factorKKT(P,A1,[],[],flag);
elseif strcmpi(flag,'ldl') % strategy similar to the direct method in SCS
    factors = factorKKT([],[],[],At,flag);
    
elseif strcmpi(flag,'blk-ldl')
    [A1,A2,D] = diviConstraint(At,K);
    factors = factorKKT(D,A1,A2,[],flag);

elseif strcmpi(flag,'blk-chol')
    factors = factorKKT([],[],[],At,flag);
end

% Compute a useful vector
z.x   = c; 
z.y   = -b;
A     = At';
eta   = solveInner(factors,At,A,z);
const = 1+c'*eta.x - b'*eta.y;
eta.x = eta.x/const; 
eta.y = eta.y/const; 

% Set function handle
solInner = @(v)solveInner(factors,At,A,v);

end

%----------------------
% divide the constraints
%----------------------
function [A1,A2,D] = diviConstraint(At,K)
    % divide the equality constraints
    % Ax = b --> [A1 A2] x = b
    % such that D = A2*A2' is diagonal
    
    A  = At';
    nCone = length(K.f) + length(K.l) + length(K.q) + length(K.s);
    nConeVars = cumsum([0,K.f, K.l, K.q, K.s.*(K.s+1)/2]);
    
    for i = nCone:-1:1     %% this strategy doesn't need to reorder the consstraints Ax = b, 
                           % but this may be improved by reordering      
        tmpA = A(:,nConeVars(i)+1:nConeVars(i+1));
        if ~isdiag(tmpA*tmpA')   %% nondiagonal part
            dFlag = i + 1;
            break
        end
        dFlag = i;
    end
    
    A1 = A(:,1:nConeVars(dFlag)); 
    A2 = A(:,nConeVars(dFlag)+1:end);
    D  = A2*A2';
    if ~isdiag(D) || isempty(D)
        error('No orthognal constraints; please use another option: hsde!')
    end 
end
 
%----------------------
% solveInner
%----------------------
function u = solveInner(factors,At,A,v)
    % Solve system of form
    %
    % [  I  -A'] [u1] = [v1]
    % [  A  I  ] [u2]   [v2]
    %
    % v.x  = v1
    % v.y  = v2
    %
    % Output is formatted the same way
    
    persistent useBuiltin
    if(isempty(useBuiltin))
        %default to look for CSparse code
         useBuiltin = ~exist(['+cdcs_utils',filesep,'cs_ltsolve.' mexext],'file');
         useBuiltin = useBuiltin | ~exist(['+cdcs_utils',filesep,'cs_lsolve.' mexext],'file');
    end
    
    % Import functions
    import cdcs_utils.cs_lsolve
    import cdcs_utils.cs_ltsolve
    
    if isfield(factors,'diagFlag') && factors.diagFlag == true   %% AA' is diagonal
        P    = factors.P;
        u.y  = P.*(v.y - A*v.x);
        u.x  = v.x + At*u.y;
    
    elseif strcmpi(factors.flag,'blk')
        %Cholesky factors and original blocks
        R   = factors.R;
        A1  = factors.A1;
        A1t = factors.A1t;
        P   = factors.P;
        s   = factors.s;
        si  = factors.si;
        
        % First find v0 to solve system like M*u2 = v0
        v2 = -(A*v.x) + v.y;
        
        % Solve system (I+AA')*u2=v2, i.e., (P+A1A1')u2 = v2 using factors
        z = P.*v2;
        d = A1t*z;
        ds= full(d(s,1));                             % permute
        if(useBuiltin) 
            p = R'\((R\ds));                    %Native matlab version (slow)
        else
            p = cs_ltsolve(R,cs_lsolve(R,ds));  %Csparse version (avoids transpose)
        end
        p = p(si,1);                            % permute back
        
        % finish solve
        u.y  = z - P.*(A1*p);
        u.x  = v.x + At*u.y;

        
    elseif strcmpi(factors.flag,'ldl')
        
        %LDL factors and permutation 
        L  = factors.L;
        D  = factors.D;
        s  = factors.s;
        si = factors.si;

        %assemble the right hand side
        b = [v.x;v.y];

        %reorder using ldl permutation
        bs = b(s,1);

       if(useBuiltin)
            %Native matlab version (slow)
            q = L'\(D\(L\bs));
       else
            %Csparse version (avoids transpose)
            q = cs_ltsolve(L,D\cs_lsolve(L,full(bs)));
       end

        %permute back
        x = q(si,1);

        %repartition into (x,y)
        nz = length(v.x);
        u.x  = x(1:nz,1);
        u.y  = -x((nz+1):end,1);
        
   elseif strcmpi(factors.flag,'blk-ldl')
        
        %LDL factors and permutation 
        m   = factors.m;
        A2  = factors.A2;
        A2t = factors.A2t;
        L   = factors.L;
        D   = factors.D;
        s   = factors.s;
        si  = factors.si;

        %assemble the right hand side 
        b  = [v.x(1:m,1);v.y-A2*v.x(m+1:end,1)];

        %reorder using ldl permutation
        bs = b(s,1);

       if(useBuiltin)
            %Native matlab version (slow)
            q = L'\(D\(L\bs));
       else
            %Csparse version (avoids transpose)
            q = cs_ltsolve(L,D\cs_lsolve(L,full(bs)));
       end

        %permute back
        x = q(si,1);

        %repartition into (x,y)
        u.x = zeros(length(v.x),1);
        u.x(1:m)  = x(1:m,1);
        
        u.y  = -x((m+1):end,1);
        u.x(m+1:end)  = A2t*u.y+v.x(m+1:end,1);
        
    elseif strcmpi(factors.flag,'blk-chol')
        
        % Cholesky factors and permutation 
        A  = factors.A;
        At = factors.At;
        R   = factors.R;
        s   = factors.s;
        si  = factors.si;

        %assemble the right hand side 
        b  = -(A*v.x) + v.y;

        %reorder using ldl permutation
        bs = b(s,1);

       if(useBuiltin)
            %Native matlab version (slow)
            q = R'\(R\bs);
       else
            q = cs_ltsolve(R,cs_lsolve(R,full(bs)));  %Csparse version (avoids transpose)
       end

        %permute back
        x = q(si,1);

        %repartition into (x,y)
        u.y = x;
        u.x = v.x + At*u.y;
    else
        error('Unknown method')
    end
    
    
end


%----------------------
% factorKKT
%----------------------
function factors = factorKKT(P,A1,A2,At,flag)

    if(nargin < 2)
        A1 = [];
        flag = 'blk';
    end

    if(nargin < 3)
        flag = 'blk';
    end
    
    if strcmpi(flag,'blk')   
        if isempty(A1)  %% no free variable,        
            factors.diagFlag = true;  %% all of the constraints are orthogonal
            factors.P        = P;     %% diagnal inverse  
        else  %% NOT all of the constraints are orthogonal, only part of them are!
            factors.diagFlag = false;
            
            [n,m] = size(A1);
            B = spdiags(sqrt(P),0,n,n)*A1;      % slow with sparse type with few nonzeros?
            M = speye(m) + sparse(B'*B);        % = I + A1'*P*A
            [R,p,s] = chol(M,'lower','vector');

            if p==0
                %For maximum efficiency in projection, store both
                %the permutation s and its inverse permutation
                factors.flag  = 'blk';
                factors.R     = R;
                factors.A1    = A1;
                factors.A1t   = A1';
                factors.P     = P;
                factors.s     = s;
                tmp           = 1:length(s);
                factors.si(s) = tmp; %inverse permutation
            end
        end
        
    elseif strcmpi(flag,'ldl')           
         % use LDL
         % Use only the lower part of the matrix for factorization
        [n,m] = size(At);
        M = [speye(n,n), sparse(At)
             sparse(m,n), -speye(m,m)];
        [U,D,s] = ldl(M,'upper','vector');
        
        %For maximum efficiency in projection, store both
        %the permutation s and its inverse permutation
        factors.flag = 'ldl';
        factors.L     = U';
        factors.D     = D;
        factors.s     = s;
        tmp           = 1:length(s);
        factors.si(s) = tmp; %inverse permutation

    elseif strcmpi(flag,'blk-ldl')           
         % first block elimination and then use LDL decomposition
         % Use only the lower part of the matrix for factorization
        [n,m] = size(A1);
        M = [speye(m,m), sparse(A1');
             sparse(n,m), -speye(n,n)-P];
        [U,D,s] = ldl(M,'upper','vector');
        
        %For maximum efficiency in projection, store both
        %the permutation s and its inverse permutation
        factors.m     = m;
        factors.A2    = A2;
        factors.A2t   = A2';
        factors.flag = 'blk-ldl';
        factors.L     = U';
        factors.D     = D;
        factors.s     = s;
        tmp           = 1:length(s);
        factors.si(s) = tmp; %inverse permutation
        
    elseif strcmpi(flag,'blk-chol')           
         % first block elimination and then use cholesky decomposition
         % Use only the lower part of the matrix for factorization
        factors.diagFlag = false;
            [~,m] = size(At);
            M = speye(m) + sparse(At'*At);        % = I + A*P*A'
            [R,p,s] = chol(M,'lower','vector');

            if p==0
                %For maximum efficiency in projection, store both
                %the permutation s and its inverse permutation
                factors.flag  = 'blk-chol';
                factors.R     = R;
                factors.A     = At';
                factors.At    = At;
                factors.P     = P;
                factors.s     = s;
                tmp           = 1:length(s);
                factors.si(s) = tmp; %inverse permutation
            end    
    else
        error('Unknown method')
    end
end
