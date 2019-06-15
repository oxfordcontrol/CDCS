function [x,y,z,info,opts] = polish(X,Y,Z,others,Kold,At,b,c,Ech,chstuff,info,opts)
    % Polish the solution from ADMM using Newtwon's steps
    % 
    %  [ -mu*Hessian   A'][delta_x] = -[C + mu*graident - A'y]
    %  [ A             0 ][delta_y] = -[      Ax - b         ]
    %
    %  2019/06/14: only aim for primal method now and single dense PSD cone
    
    % test 
    
    gH_Params.n   = Kold.s*(Kold.s+1)/2; 
    gH_Params.Q   = svecTransMat(Kold.s); 
    gH_Params.nu  = Kold.s;
    
    xInit = svec(X.blk{1}+1e-3*eye(Kold.s)); % make sure the solution is PSD. 
    yInit = others.dual;
    
    mu = opts.imp.tol/gH_Params.nu;  % this is not good
    
    xCurrent = xInit; yCurrent = yInit;
    [n,m] = size(At);
    for iter = 1:opts.ipm.maxIter
        
        % Form Newton system
        [in, g, H, L] = gH_psd(xCurrent, gH_Params);
        RHS = [-mu*H  At;
                At'  zeros(m,m)];
                  
        LHS = [c + mu*g - At*yCurrent;b - At'*xCurrent];
        
        % linear solver
        delta = RHS\LHS;
        
        % Compute step length via line search
        
        
        % Update
        
        % other stuff
               
    end
    
end

