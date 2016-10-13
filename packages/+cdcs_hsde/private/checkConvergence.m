function [Flag,pcost,dcost,presi,dresi,gap,opts] = checkConvergence(u,v,iter,admmtime,opts,At,b,c,E,hatu,uold)

% CHECKCONVERGENCE primal/dual infeasible

Flag.stop = false;
Flag.isConverged = false;   %% converged
Flag.pinfeas     = false;   %% primal infeasible
Flag.dinfeas     = false;   %% dual infeasible

if u.tau > 0   %% feasible problem
    x = u.x./u.tau; y = u.y./u.tau;
    dcost = (b'*y)/opts.scaleFactors.sc_cost;        %% dual cost
    pcost = (c'*x)/opts.scaleFactors.sc_cost;        %% primal cost
    
    % this is based on the scaled data
    presi = norm(At'*x - b,'fro')./max([1,opts.nAt,opts.nb]);
    %dresi = norm((At*y + accumarray(E,u.v./u.tau.*opts.scaleFactors.E2).*opts.scaleFactors.D1- c),'fro')./max([1,opts.nAt,opts.nc]);
    dresi = norm((At*y + accumarray(E,u.v./u.tau)- c),'fro')./max([1,opts.nAt,opts.nc]);
    gap   = abs(pcost - dcost)/(1 + abs(pcost) + abs(dcost));

    % this is based on the original data, see the paper by O'Donoghue et al.
%     presi = norm((At'*x - b)./opts.scaleFactors.E1/opts.scaleFactors.sc_b,'fro')./(1+opts.nub);
%     dresi = norm((At*y + accumarray(E,u.v./u.tau.*opts.scaleFactors.E2).*opts.scaleFactors.D1- c)./opts.scaleFactors.D1/opts.scaleFactors.sc_c,'fro');
%     dresi = dresi./(1+opts.nuc);
%     gap   = abs(pcost - dcost)/(1 + abs(pcost) + abs(dcost));
    
    if max([presi,dresi,gap]) < opts.relTol
        Flag.isConverged = true;
        Flag.stop = true;
    end
    
    %% adpative penalty?
    % GF: Only adapt penalty if not stopping after this iteration...
    if opts.adaptive && ~Flag.stop
        
        % standard ADMM residual according to Boyd's paper
        utmp = vertcat(u.x,u.xh,u.y,u.v,u.tau); 
        vtmp = vertcat(v.x,v.xh,v.y,v.v,v.kappa);
        hatutmp = vertcat(hatu.x,hatu.xh,hatu.y,hatu.v,hatu.tau);
        uoldtmp = vertcat(uold.x,uold.xh,uold.y,uold.v,uold.tau);
        
        r = norm(utmp-hatutmp,'fro')./max(norm(utmp,'fro'),norm(hatutmp,'fro'));
        s = opts.rho*norm(utmp-uoldtmp,'fro')./max(norm(utmp,'fro'),norm(vtmp,'fro'));
        
        % Update ADMM penalty parameter
        opts = updatePenalty(opts,r,s);
        
    end
    
else % infeasible or unbounded problem?
    pinfIndex = b'*u.y;  % index of certificating primal infesibility
    dinfIndex = -c'*u.x; % index of certificating dual infesibility
    if dinfIndex > 0    % dual infeasible
        pfeas    = norm(At'*u.x,'fro');
        pfeasTol = dinfIndex/opts.nc*opts.relTol;
        if pfeas <= pfeasTol  %% a point x that certificates dual infeasiblity
           Flag.dinfeas = true;
           Flag.stop = true;
        end
    end
    if pinfIndex > 0   % primal infeasible
        dfeas    = norm(At*u.y+accumarray(E,u.v),'fro');
        dfeasTol = pinfIndex/opts.nb*opts.relTol;
        if dfeas <= dfeasTol  %% a point y that certificates primal infeasiblity
           Flag.pinfeas = true;
           Flag.stop = true;
        end
    end
        
    %% value
    presi = NaN;
    dresi = NaN;
    pcost = NaN;
    dcost = NaN;
    gap   = NaN;
end

%progress message
if opts.verbose && (iter == 1 || ~mod(iter,opts.dispIter) || Flag.stop)
    fprintf('%5d | %7.2e | %7.2e | %9.2e | %9.2e | %8.2e | %8.2e | %8.2e |\n',...
            iter,presi,dresi,pcost,dcost,gap, opts.rho,toc(admmtime))
end  

end

%% Nested functions

function opts = updatePenalty(opts,pres,dres)
% Update penaly

    persistent itPinf itDinf
    if isempty(itPinf) || isempty(itDinf)
        % initialize persistent iteration counters when entering for the first
        % time (persistent variables are empty and not zero when declared)
        itPinf = 0; % number of iterations for which pinf/dinf <= eta
        itDinf = 0; % number of iterations for which pinf/dinf > eta
    end

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

end

