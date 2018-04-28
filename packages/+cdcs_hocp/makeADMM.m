function [step1,step2,step3,checkConv] = makeADMM(At,b,c,K,Ech,opts)

% Make ADMM operators

% divide equally
IA  = accumarray(Ech,1); % d
cd  = c./IA; cd  = cd(Ech);

% Useful stuff
[projAffine,projCone] = makeProjectors(At,b,c,K,cd,Ech,opts);
                        % makeProjectors(At,b,c,K,cd,E,opts)

% Update steps
step1 = @(X,Y,Z,rho,others)updateX(X,Y,Z,rho,others,Ech,K,projAffine,opts);
step2 = @(X,Y,Z,rho,others)updateY(X,Y,Z,rho,others,projCone);
step3 = @(X,Y,Z,rho,others)updateZ(X,Y,Z,rho,others,K);

% Convergence check
checkConv = @(X,Y,Z,YOld,others,iter,admmtime,opts,log)...
    checkConvergence(X,Y,Z,YOld,others,b,c,Ech,iter,opts,admmtime,log);



end


