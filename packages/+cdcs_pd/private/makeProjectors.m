function [projAffine,projCone] = makeProjectors(At,b,c,K,cd,E,opts)
% Create projection operators

% Import functions
import cdcs_utils.makeConeVariables
import cdcs_utils.blockify

%rho = opts.rho;
H = accumarray(E,1);

% split cost between sub-variables?
%----------------------------------
if strcmpi(opts.solver,'dual') || ~opts.yPenalty
    cd  = sparse(opts.nXk,1); %all zeros unless primal form & zPenalty options
end
CD = makeConeVariables(K);
CD = blockify(CD,cd,K);

% Projection on affine subspace
%------------------------------
factors = factorKKT(H,At,full(b),opts.KKTfact);
projAffine = @(z,lambda,rho)projectAffine(factors,z,b,c,lambda,E,rho,opts);

% The projection onto K
%----------------------
projCone  = @(EX,L,rho)projectionToCones(K,EX,L,CD,rho,opts);


% END FUNCTION
end


%----------------------
% projectAffine
%----------------------
function [x,y] = projectAffine(factors,z,b,c,lambda,E,rho,opts)
    
    % Find the optimal (x,y) pair by solving a KKT system
    if strcmpi(opts.solver,'primal')
        rx = full( accumarray(E,lambda./rho + z) - c./rho );    %compute the RHS
        [x,y] = solveKKT(factors,rx,full(b));                   %solve the KKT system 
        if opts.yPenalty
            % account for factor of 2 due to penalty
            y = -y.*(rho/2);        
        else
            y = -y.*rho;
        end
        
    elseif strcmpi(opts.solver,'dual')
        rx = full( c - accumarray(E,lambda./rho + z) );    %compute the RHS
        [x,y] = solveKKT(factors,rx,-full(b)./rho);                   %solve the KKT system 
        x = -x.*rho;
        
    end

end


%----------------------
% projectionToCones
%----------------------
function S = projectionToCones(K,EX,L,CZ,rho,opts)
    
    % Import functions
    import cdcs_utils.projectK

    % project variables to cones
    if strcmpi(opts.solver,'primal')
        S = cellfun(@(EX,L,CZ)(EX - L./rho - CZ./rho),EX,L,CZ,'UniformOutput',false);
        S = projectK(S,K,0);
        
    elseif strcmpi(opts.solver,'dual')
        S = cellfun(@(EX,L)(EX - L./rho),EX,L,'UniformOutput',false);
        S = projectK(S,K,1);
        
    end

   

end



%----------------------
% factorKKT
%----------------------
function factors = factorKKT(H,At,b,flag)

    % Compute factors for solving the KKT system
    %
    %  [ Diag(H)  A'] [z] = [bz]
    %  [      A   0 ] [p]   [bp]
    %
    % where H is a vector. Direct factorization if flag='ldl', otherwise use block
    % elimination or inversion.

    if(nargin < 2)
        At = [];
        flag = 'blk';
    end

    if(nargin < 3)
        flag = 'blk';
    end
    
    [n,m] = size(At);

    if strcmpi(flag,'ldl')

        %Use only the lower part of the matrix for factorization
        M = [spdiags(H,0,n,n), sparse(At)
            sparse(m,n), sparse(m,m)];
        [U,D,s] = ldl(M,'upper','vector');


        %For maximum efficiency in projection, store both
        %the permutation s and its inverse permutation
        factors.flag = 'ldl';
        factors.L     = U';
        factors.D     = D;
        factors.s     = s;
        tmp           = 1:length(s);
        factors.si(s) = tmp; %inverse permutation


    elseif strcmpi(flag,'inv')

        % Invert matrix!
        B = spdiags(1./sqrt(H),0,n,n)*At;   % slow with sparse type with few nonzeros?
        Q = sparse(B'*B);                   % = A*inv(H)*A'
        Q = Q\speye(m);

        %For max speed, check if Q is sparse or not
        if nnz(Q)/(m^2) > 0.6;
            Q = full(Q);
        else
            Q = sparse(Q);
        end
        factors.flag  = 'inv';
        factors.QA    = Q*(At');
        factors.Qbp   = Q*b;
        factors.At    = At;
        factors.H     = H;


    else

        B = spdiags(1./sqrt(H),0,n,n)*At;   % slow with sparse type with few nonzeros?
        M = sparse(B'*B);                   % = A*inv(H)*A'
        [R,p,s] = chol(M,'lower','vector');

        if p==0
            %For maximum efficiency in projection, store both
            %the permutation s and its inverse permutation
            factors.flag = 'blk';
            factors.R     = R;
            factors.A     = At';
            factors.At    = At;
            factors.H     = H;
            factors.s     = s;
            tmp           = 1:length(s);
            factors.si(s) = tmp; %inverse permutation
        else
            % not definite - use LDL
            M = [spdiags(H,0,n,n), sparse(At)
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


%----------------------
% solveKKT
%----------------------
function [z,p] = solveKKT(factors,bz,bp)

    % Find a solution for the KKT system
    %
    %  [ diag(H)  A'] [z] = [bz]
    %  [      A   0 ] [p]   [bp]
    %
    %  given an appropriate factor structure

    %NB : The below is very inefficient in the case
    %that the rank-1 perturbations is baked in.  
    %better to do the solve with the factors of
    %H and then apply rank-1 updates for uu'

    %Two strategies:  If H is diagonal, then do block reduction.
    %otherwise, factor the complete LHS above

    persistent useBuiltin
    if(isempty(useBuiltin))
        %default to look for CSparse code
         useBuiltin = ~exist(['+cdcs_utils',filesep,'cs_ltsolve.' mexext],'file');
         useBuiltin = useBuiltin | ~exist(['+cdcs_utils',filesep,'cs_lsolve.' mexext],'file');
    end
    
    % Import functions
    import cdcs_utils.cs_lsolve
    import cdcs_utils.cs_ltsolve
    
    % Compute

    if strcmpi(factors.flag,'ldl')

        %LDL factors and permutation 
        L  = factors.L;
        D  = factors.D;
        s  = factors.s;
        si = factors.si;

        %assemble the right hand side
        b = [bz;bp];

        %reorder using Cholesky permutation
        bs = b(s,1);

       if(useBuiltin)
            %Native matlab version (slow)
            q = L'\(D\(L\bs));
       else
            %Csparse version (avoids transpose)
            q = cs_ltsolve(L,D\cs_lsolve(L,bs));
       end

        %permute back
        x = q(si,1);

        %repartition into (z,p)
        nz = length(bz);
        z  = x(1:nz,1);
        p  = x((nz+1):end,1);


    elseif strcmpi(factors.flag,'inv')

        % Inverses and product matrices
        H = factors.H;
        QA = factors.QA;
        Qbp = factors.Qbp;
        At = factors.At;

        p = QA*(bz./H) - Qbp;
        z = (bz - At*p)./H;

    elseif strcmpi(factors.flag,'blk')

        %Cholesky factors and original blocks
        R = factors.R;
        A = factors.A;
        At = factors.At;
        H = factors.H;

        s  = factors.s;
        si = factors.si;

        w = bz./H;
        d = (A*w) - bp;

        %reorder using Cholesky permutation
        ds = d(s,1);

        if(useBuiltin)
            %Native matlab version (slow)
            p = R'\((R\ds));
        else
            %Csparse version (avoids transpose)
            p = cs_ltsolve(R,cs_lsolve(R,ds));
        end

        %permute back
        p = p(si,1);

        % Finish solve
        z = w - (At*p)./H;

    end

end
    
    