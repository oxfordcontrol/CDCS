function [xi,solInner] = factorMatrix_scale(At,b,c,E,flag,scaleFactors)

% generate projector onto affine constraints
% H = E2*H*D1, F = D2*E2

scaleD1 = scaleFactors.D1;   
scaleD2 = scaleFactors.D2;
scaleE2 = scaleFactors.E2;

F = scaleD2.*scaleE2;

D = accumarray(E,scaleE2.^2.*(1-F./(1+F.^2).*F)).*(scaleD1.^2);
P = 1./(1+D);

% First factor the matrix (I+APA') where P = inv(I+D) diagonal
factors = factorKKT(P,At,flag);

% Compute a useful vector
z.x  = c; 
z.xh = sparse(length(E),1);
z.y  = -b;
z.v  = sparse(length(E),1);
xi = solveInner(factors,E,z,scaleD1,scaleE2,F);
const = 1+c'*xi.x - b'*xi.y;
xi.x  = xi.x/const; 
xi.xh = xi.xh/const; 
xi.y  = xi.y/const; 
xi.v  = xi.v/const;

% Set function handle
solInner = @(v)solveInner(factors,E,v,scaleD1,scaleE2,F);

end

%----------------------
% solveInner
%----------------------
function u = solveInner(factors,E,v,scaleD1,scaleE2,scaleF)
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
    
    if strcmpi(factors.flag,'blk')

        invF = 1./(1+scaleF.^2);
        %Cholesky factors and original blocks
        R  = factors.R;
        A  = factors.A;
        At = factors.At;
        P  = factors.P;
        s  = factors.s;
        si = factors.si;
        
        % First find v0 to solve system like M*u2 = v0
        v01 = v.x + accumarray(E,v.v.*scaleE2).*scaleD1 + (At*v.y) ;
        v02 = v.xh - v.v.*scaleF;
        
        % Solve system M*u2=v0 using factors
        z = P.*(v01 + accumarray(E,scaleE2.*scaleF.*invF.*v02).*scaleD1);
        d = A*z;
        ds= d(s,1);                             % permute
        if(useBuiltin) 
            p = R'\((R\ds));                    %Native matlab version (slow)
        else
            p = cs_ltsolve(R,cs_lsolve(R,full(ds)));  %Csparse version (avoids transpose)
        end
        p = p(si,1);                            % permute back
        
        % finish solve
        u.x  = z - P.*(At*p);
        tmp = u.x.*scaleD1;
        u.xh = v02.*invF + invF.*scaleF.*scaleE2.*tmp(E);
        u.y  = v.y - A*u.x;
        u.v  = v.v + u.xh.*scaleF - scaleE2.*tmp(E);
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
    
    [n,m] = size(At);

    if strcmpi(flag,'ldl')
        % TO DO


    elseif strcmpi(flag,'inv')
        % TO DO
        


    else
               
        [~,m] = size(At);
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
            % not definite - use LDL
            M = [spdiags(1./P,0,n,n), sparse(At)
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
