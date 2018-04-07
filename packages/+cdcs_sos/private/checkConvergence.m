function [stop,info,log,opts] = checkConvergence(hatu,u,v,uold,iter,admmtime,opts,At,b,c,others,log)
% CHECKCONVERGENCE primal/dual infeasible
%
% Problem codes
% 0: converged
% 1: primal infeasible
% 2: dual infeasible
% 3: max number of iterations reached

stop         = false;
info.problem = 3;

if u.tau > 0   %% feasible problem
    x = u.x./u.tau;
    y = u.y./u.tau;
    dcost = (b.'*y)/opts.scaleFactors.sc_cost;        %% dual cost
    pcost = (c.'*x)/opts.scaleFactors.sc_cost;        %% primal cost
      
    % Residuals before rescaling (O'Donoghue et al, Section 5)
    % Use same normalization as O'Donoghue et al
    pk  = (At.'*x - b).*opts.scaleFactors.E;
    presi = ( norm(pk,'fro')./(1+opts.nb_init) ) / opts.scaleFactors.sc_b;
    dk  = (At*y + (v.x./u.tau) - c).*opts.scaleFactors.D;
    dresi = ( norm(dk,'fro')./(1+opts.nc_init) ) / opts.scaleFactors.sc_c;
    gap = abs(pcost - dcost)/(1 + abs(pcost) + abs(dcost));
   
    if max([presi,dresi,gap]) < opts.relTol
        info.problem = 0;
        stop         = true;
    end
   
else % infeasible or unbounded problem?
    pinfIndex = b.'*u.y;  % index of certificating primal infesibility
    dinfIndex = -c.'*u.x; % index of certificating dual infesibility
    if dinfIndex > 0    % dual infeasible
        % Scaled variables
%         pfeas    = norm(At'*u.x,'fro');
%         pfeasTol = dinfIndex/opts.nc*opts.relTol;
        % Original variables
        pfeas    = norm( (At'*u.x).*opts.scaleFactors.E ,'fro');
        pfeasTol = dinfIndex/opts.nc_init/opts.scaleFactors.sc_c*opts.relTol;
        if pfeas <= pfeasTol  %% a point x that certificates dual infeasiblity
            info.problem = 2;
            stop = true;
        end
    end
    if pinfIndex > 0  % primal infeasible
        % Original variables
        dfeas    = norm( (At*u.y+v.x).*opts.scaleFactors.D ,'fro');
        dfeasTol = pinfIndex/opts.nb_init/opts.scaleFactors.sc_b*opts.relTol;
        if dfeas <= dfeasTol  %% a point y that certificates primal infeasiblity
            info.problem = 1;
            stop = true;
        end
    end
    
    % value
    presi = NaN;
    dresi = NaN;
    pcost = NaN;
    dcost = NaN;
    gap   = NaN;
end

%progress message
if opts.verbose && (iter == 1 || ~mod(iter,opts.dispIter) || stop)
    if ~isnan(presi)
        fprintf('%5d | %7.2e | %7.2e | %9.2e | %9.2e | %8.2e |  %8.2e |\n',...
            iter,presi,dresi,pcost,dcost,gap, toc(admmtime))
    else
        % for the purpose of alginment
        fprintf('%5d | %8d | %8d | %9d | %9d | %8d |  %8.2e |\n',...
            iter,presi,dresi,pcost,dcost,gap, toc(admmtime))
    end
        
end

% log information
% Use preallocation for speed
if iter==1
	% Much faster preallocation method
    log.dcost = zeros(opts.maxIter,1);
    log.gap   = zeros(opts.maxIter,1);
end

log.pres(iter)  = presi;
log.dres(iter)  = dresi;
log.cost(iter)  = pcost;   %% use primal cost
log.dcost(iter) = dcost;
log.gap(iter)   = gap;

% end main
end

