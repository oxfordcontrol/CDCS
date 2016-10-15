function [stop,info,log,opts] = checkConvergence(X,Y,Z,YOld,others,b,c,E,iter,opts,admmtime)

% CHECKCONVERGENCE
% Use the basic convergence test in the Boyd survey paper

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
    info.problem = 1;
end

%progress message
if opts.verbose && (iter == 1 || ~mod(iter,opts.dispIter) || stop)
    fprintf('%5d | %8.2e | %8.2e | %9.2e  | %8.2e | %8.2e |\n',...
        iter,pres,dres,cost,opts.rho,toc(admmtime));
end


% Update penalty parameter
if opts.adaptive
    resRat = pres/dres;
    if resRat >= opts.mu
        itPinf = itPinf+1;
        itDinf = 0;
        if itPinf >= opts.rhoIt
            % ratio of pinf and dinf remained large for long => rescale rho
            itPinf = 0;
            opts.rho = min(opts.rho*opts.tau, opts.rhoMax);
        end
    elseif 1/resRat >= opts.mu
        itDinf = itDinf+1;
        itPinf = 0;
        if itDinf >= opts.rhoIt
            % ratio of pinf and dinf remained small for long => rescale rho
            itDinf = 0;
            opts.rho = max(opts.rho/opts.tau, opts.rhoMin);
        end
    end
end

% Log errors
log(iter).pres = pres;
log(iter).dres = dres;
log(iter).cost = cost;

end