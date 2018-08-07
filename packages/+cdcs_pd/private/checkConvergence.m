function [stop,info,log,opts] = checkConvergence(X,Y,Z,YOld,others,b,c,E,iter,opts,admmtime,log)

% CHECKCONVERGENCE
% Use the basic convergence test in the Boyd survey paper
% Primal/dual solver cannot detect infeasibility so error codes can only be
% 0: problem successfully solved
% 3: max number of iterations reached


persistent itPinf itDinf
if iter == 1
    % initialize persistent iteration counters when entering for the first
    % time (persistent variables are empty and not zero when declared)
    itPinf = 0; % number of iterations for which pinf/dinf <= eta
    itDinf = 0; % number of iterations for which pinf/dinf > eta
end

% Import functions
import cdcs_utils.flatten

% Extract some variables
Ex = flatten(Y.vec,X.blk);
y = Y.vec;
yOld = YOld.vec;
z = Z.vec;
rho = opts.rho;

% primal residual
r = norm(Ex-y,'fro');
pres  = r./max(norm(Ex,'fro'),norm(y,'fro'));

if strcmpi(opts.solver,'primal')
    % cost
    x = X.vec;
    cost = full(c.'*x)/opts.scaleFactors.sc_cost;
    
    % dual residual
    s = rho.*(norm( accumarray(E,(y-yOld)) ,'fro'));
    dres  = s./norm(accumarray(E,z),'fro');
    
elseif strcmpi(opts.solver,'dual')
    % cost
    x = X.vec(1:opts.m);
    cost = full(b.'*x)/opts.scaleFactors.sc_cost;
    
    % dual residual
    s = rho.*(norm(y-yOld,'fro'));
    dres  = s./norm(z,'fro');
end

% Costs


%stopping criteria
if(max(pres,dres)<opts.relTol)
    stop = true;
    info.problem = 0;
else
    stop = false;
    info.problem = 3;
end

%progress message
if opts.verbose && (iter == 1 || ~mod(iter,opts.dispIter) || stop)
    fprintf('%5d | %8.2e | %8.2e | %9.2e  | %8.2e | %8.2e |\n',...
        iter,pres,dres,cost,opts.rho,toc(admmtime));
end


% Update penalty parameter
if opts.adaptive
    resRat = pres/dres;
    if resRat >= opts.nu
        itPinf = itPinf+1;
        itDinf = 0;
        if itPinf >= opts.rhoIt
            % ratio of pinf and dinf remained large for long => rescale rho
            itPinf = 0;
            opts.rho = min(opts.rho*opts.mu, opts.rhoMax);
        end
    elseif 1/resRat >= opts.nu
        itDinf = itDinf+1;
        itPinf = 0;
        if itDinf >= opts.rhoIt
            % ratio of pinf and dinf remained small for long => rescale rho
            itDinf = 0;
            opts.rho = max(opts.rho/opts.mu, opts.rhoMin);
        end
    end
end

% Log errors
log.pres(iter) = pres;
log.dres(iter) = dres;
log.cost(iter) = cost;

%=======
% Use preallocation for speed
%if iter==1
%    cc = cell(opts.maxIter,1);
%    log = struct('pres',cc,'dres',cc,'cost',cc);
%end
%log(iter).pres = pres;
%log(iter).dres = dres;
%log(iter).cost = cost;
%>>>>>>> master

end