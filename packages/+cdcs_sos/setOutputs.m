function [x,y,z,info,opts] = setOutputs(X,Y,Z,K,info,opts)

% SETOUTPUTS.M
%
% Set outputs in sedumi format using positive semidefinite completion algorithms

% Import functions
import cdcs_utils.makeConeVariables
import cdcs_utils.flatten
import cdcs_utils.blockify

x = zeros(opts.n_init,1);
y = zeros(opts.m_init,1);
z = zeros(opts.n_init,1);


%% vector value
xtmp = Y.x/Y.tau;
ytmp = Y.y/Y.tau;
ztmp = Z.x/X.tau; 

% Scale solution
% Scale before blockify because rescaling was done after svec operation in
% rescaleData.m
xtmp = (xtmp./opts.scaleFactors.D)./opts.scaleFactors.sc_b;
ytmp = (ytmp./opts.scaleFactors.E)./opts.scaleFactors.sc_c;
ztmp = (ztmp.*opts.scaleFactors.D)./opts.scaleFactors.sc_c;



% Reorder the variables x, z back to match the original data 
% becasue we reordered the data into the prepcosseing step
Kreorder.s = K.s(opts.sos.ReOrder);                        % This is the PSD cone after reordering
nConeVars = cumsum([K.f+K.l+K.q, Kreorder.s.*(Kreorder.s+1)/2]);
[~,ReOrder] = sort(opts.sos.ReOrder);               %% reorder the variables x
Index = 1;
ReOrderVariables = zeros(nConeVars(end)-nConeVars(1),1);
for i = 1:length(K.s)
    Num = nConeVars(ReOrder(i)+1) - nConeVars(ReOrder(i));
    ReOrderVariables(Index:Index+Num-1)  = nConeVars(ReOrder(i))+1:nConeVars(ReOrder(i)+1);
    Index = Index + Num;
end
ReOrderVariables = [(1:nConeVars(1))';ReOrderVariables];
xtmp = xtmp(ReOrderVariables);
ztmp = ztmp(ReOrderVariables);


% block value
[xmat,zmat] = makeConeVariables(K);        % blockified
xmat  = blockify(xmat,xtmp,K);
zmat  = blockify(zmat,ztmp,K);
xtmp = flatten(xmat,xmat,0); 
ztmp = flatten(zmat,zmat,0); 

% % Scale solution
% xtmp = (xtmp./opts.scaleFactors.D)./opts.scaleFactors.sc_b;
% ytmp = (ytmp./opts.scaleFactors.E)./opts.scaleFactors.sc_c;
% ztmp = (ztmp.*opts.scaleFactors.D)./opts.scaleFactors.sc_c;

% Assign solution to used variables
x(opts.usedvars) = xtmp;
z(opts.usedvars) = ztmp;
y = ytmp;

% END FUNCTION
end