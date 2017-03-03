function [xi,solInner] = factorMatrix(At,b,c,E,flag)
% generate projector onto affine constraints
% Here, we do not consider the modification of H in the scaling process

D = accumarray(E,1);
P = 1./(1+D./2);

% First factor the matrix (I+APA') where P = inv(I+1/2D) diagonal
factors = factorKKT(P,At,flag);

% Compute a useful vector
z.x  = c; 
z.xh = sparse(length(E),1);
z.y  = -b;
z.v  = sparse(length(E),1);
xi   = solveInner(factors,E,z);
const = 1+c'*xi.x - b'*xi.y;
xi.x  = xi.x/const; 
xi.xh = xi.xh/const; 
xi.y  = xi.y/const; 
xi.v  = xi.v/const;

% Set function handle
solInner = @(v)solveInner(factors,E,v);

end

%----------------------
% solveInner
%----------------------
function u = solveInner(factors,E,v)
    % Solve system of form
    %
    % [  I  M1] [u1] = [v1]
    % [-M1'  I] [u2]   [v2]
    %
    % where u1 = [u11,u12], u2 = [u21,u22] etc. Input v is formatted so that 
    % v.x  = v11
    % v.xh = v12
    % v.y  = v21
    % v.v  = v22
    %
    % Output is formatted the same way
    
    persistent useBuiltin

    if(isempty(useBuiltin))
        %default to look for CSparse code
         useBuiltin = ~exist(['cs_ltsolve.' mexext],'file');  
    end
    
    % Import functions
    import cdcs_utils.cs_lsolve
    import cdcs_utils.cs_ltsolve
    
    if strcmpi(factors.flag,'blk')
        
        % Cholesky factors and original blocks
        R  = factors.R;
        A  = factors.A;
        At = factors.At;
        P  = factors.P;
        s  = factors.s;
        si = factors.si;
        
        % First find v0 to solve system like M*u2 = v0
        v01 = accumarray(E,v.v) + (At*v.y) + v.x;
        v02 = v.xh - v.v;
        
        % Solve system M*u2=v0 using factors
        z = P.*(v01 + accumarray(E,v02)./2);
        d = A*z;
        ds= full(d(s,1));                       % permute
        if(useBuiltin) 
            p = R'\((R\ds));                    % Native matlab version (slow)
        else
            p = cs_ltsolve(R,cs_lsolve(R,ds));  % Csparse version (avoids transpose)
        end
        p = p(si,1);                            % permute back
        
        % finish solve
        u.x  = z - P.*(At*p);
        u.xh = (v02 + u.x(E))./2;
        u.y  = v.y - A*u.x;
        u.v  = v.v + u.xh - u.x(E);
    else
        error('Unknown method')
    end
    
    
end


%----------------------
% factorKKT
%----------------------
function factors = factorKKT(P,At,flag)

    if(nargin < 2)
        At = [];
        flag = 'blk';
    end

    if(nargin < 3)
        flag = 'blk';
    end
    
    %[n,~] = size(At);

    if strcmpi(flag,'ldl')
        % TO DO
        % I don't think ldl decomposition is a good choice for hsde.
        error('Unknown method')

    elseif strcmpi(flag,'inv')
        % TO DO
        error('Unknown method')
    else
        [n,m] = size(At);
        B = spdiags(sqrt(P),0,n,n)*At;      % slow with sparse type with few nonzeros?
        M = speye(m) + sparse(B'*B);        % = I + A*P*A'
        [R,p,s] = chol(M,'lower','vector');

        if p==0
            %For maximum efficiency in projection, store both
            %the permutation s and its inverse permutation
            factors.flag = 'blk';
            factors.R     = R;
            factors.A     = At';
            factors.At    = At;
            factors.P     = P;
            factors.s     = s;
            tmp           = 1:length(s);
            factors.si(s) = tmp; %inverse permutation
        else
            error('Cholesky decomposition fails ... ')
            % not definite - use LDL
            % M = [spdiags(1./P,0,n,n), sparse(At)
            %     sparse(m,n), sparse(m,m)];
            % [U,D,s] = ldl(M,'upper','vector');
            % factors.flag = 'ldl';
            % factors.L     = U';
            % factors.D     = D;
            % factors.s     = s;
            % tmp           = 1:length(s);
            % factors.si(s) = tmp; %inverse permutation
        end
    end

end
