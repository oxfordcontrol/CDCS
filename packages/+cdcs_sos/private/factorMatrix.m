function [eta,solInner] = factorMatrix(At,b,c,K,flag)
% generate projector onto affine constraints
% Here, we do not consider the modification of H in the scaling process

% check orthogonality A = [A1 A2], A2 correponds to PSD cone variabels
A  = At';
A1 = A(:,1: K.f+K.l); 
A2 = A(:,K.f+K.l+1:end);
D  = A2*A2';

if ~isdiag(D)
    error('Constraints are not orthogonal!')
end

Ddiag = diag(D);   %% store the diagonal elements
P = 1./(1+Ddiag);

% First factor the matrix (I+A1'PA1) where P = inv(I+D) is a diagonal
% matrix; if A1 is empty (all the constraints are orthogonal), then return
% P
factors = factorKKT(P,A1,flag);

% Compute a useful vector
z.x  = c; 
z.y  = -b;
eta = solveInner(factors,At,z);
const = 1+c'*eta.x - b'*eta.y;
eta.x  = eta.x/const; 
eta.y  = eta.y/const; 

% Set function handle
solInner = @(v)solveInner(factors,At,v);

end

%----------------------
% solveInner
%----------------------
function u = solveInner(factors,At,v)
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
    % Import functions
    import cdcs_utils.cs_lsolve
    import cdcs_utils.cs_ltsolve

    if(isempty(useBuiltin))
        %default to look for CSparse code
         useBuiltin = ~exist(['cs_ltsolve.' mexext],'file');  
    end
    
    if factors.diagFlag == true   %% AA' is diagonal
        P    = factors.P;
        u.y  = P.*(v.y - At'*v.x);
        u.x  = v.x + At*u.y;
    
    elseif strcmpi(factors.flag,'blk')
        %Cholesky factors and original blocks
        R  = factors.R;
        A1  = factors.A1;
        A1t = factors.A1t;
        P  = factors.P;
        s  = factors.s;
        si = factors.si;
        
        % First find v0 to solve system like M*u2 = v0
        v2 = -(At.'*v.x) + v.y;
        
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
    else
        error('Unknown method')
    end
    
    
end


%----------------------
% factorKKT
%----------------------
function factors = factorKKT(P,A1,flag)

    if(nargin < 2)
        A1 = [];
        flag = 'blk';
    end

    if(nargin < 3)
        flag = 'blk';
    end
    
    if isempty(A1)  %% no free variable,        
        factors.diagFlag = true;  %% all of the constraints are orthogonal
        factors.P        = P;     %% diagnal inverse     
    else            %% NOT all of the constraints are orthogonal, only part of them are!
        factors.diagFlag = false;
        if strcmpi(flag,'ldl')
            % TO DO
            error('Unknown method')

        elseif strcmpi(flag,'inv')
            % TO DO
            error('Unknown method')
        else
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
            else
                % not definite - use LDL
                M = [spdiags(1./P,0,n,n), sparse(A1)
                     sparse(m,n), sparse(m,m)];
                [U,D,s] = ldl(M,'upper','vector');
                factors.flag = 'ldl';
                factors.L     = U';
                factors.D     = D;
                factors.s     = s;
                tmp           = 1:length(s);
                factors.si(s) = tmp; %inverse permutation
            end
        end
    end

end
