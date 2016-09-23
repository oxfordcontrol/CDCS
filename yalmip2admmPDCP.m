function [At,b,c,K] = yalmip2admmPDCP(F,obj,dual)

% YALMIP2ADMMPDCP
%
% [At,b,c,K] = yalmip2admmPDCP(F,obj) creates inputs for admmPDCP from a YALMIP 
%              conic optimization problem with objective 'obj' and constraint 
%              structure 'F' in a suitable format. The primal or dual problem is
%              exported from YALMIP according to the optional input flag 
%              'dualize'.
%
%               dualize = 0 (default): export YALMIP's primal form
%               dualize = 1          : export YALMIP's dual form
%
% Example. Create a feasible random SDP with one equality constraint in the form
%          
%                          min   trace(C*X)
%           (1)     subject to   trace(A*X) = 0,
%                                X >= 0.
%
%          in YALMIP. First, we create a primal dual solution:
%
%          >> Xopt = speye(25); yopt = 3; Zopt = 2*speye(25);
%       
%          Then we set up the data to get a feasible SDP with the solution we
%          just created:
%
%          >> A = sprandsym(25,0.35);           % a random data matrix
%          >> b = trace(A*Xopt);                % feasible primal constraint
%          >> C = A*yopt + Zopt;                % feasible dual constraint
%
%          Finally,  let us set up the optimization problem (1) in YALMIP:
%
%          >> X = sdpvar(25);                   % the variable
%          >> F = [trace(A*X)==b; X>=0];        % the constraints
%          >> obj = trace(C*X);                 % the objective   
%
%           Since problem (1) above is in YALMIP's standard dual form, we export
%           it to admmPDCP and solve it with
%
%           >> [At,b,c,K] = yalmip2admmPDCP(F,obj,1);
%           >> admmPDCP(At,b,c,K)
%
%
% See also admmPDCP

% Dualize 
if nargin==3 && dual
    [F,obj] = dualize(F,obj);
end

% Export using sedumi format
mod = export(F,obj,sdpsettings('solver','sedumi'));

% Extract data matrices and supported cones
At = mod.A;
b = mod.b;
c = mod.C;
K.f = mod.K.f;
K.l = mod.K.l;
K.q = mod.K.q;
K.s = mod.K.s;

% Check if other cones at all - cannot solve then!
if ~all(mod.K.c)==0 || ~all(mod.K.r)==0 || ~all(mod.K.p)==0
    error('Unsupported cone constraint types.')
end
