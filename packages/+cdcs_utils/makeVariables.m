function [X,Y,Z,others] = makeVariables(K,initVars,opts)

% CDCS/packages/+cdcs_utils/MAKEVARIABLES.m
%
% Initialize variables for CDCS. Call the correct variable initialization
% routine according to the method specified in opts.solver.

switch lower(opts.solver)
   
    case {'primal', 'dual'}
        [X,Y,Z,others] = cdcs_pd.makeVariables(K,initVars,opts);
        
    case {'hsde'}
        [X,Y,Z,others] = cdcs_hsde.makeVariables(K,initVars,opts);
        
    case {'sos'}
        [X,Y,Z,others] = cdcs_sos.makeVariables(K,initVars,opts);    
        
    otherwise
        error('Unknown value for ''options.solver''.')
            
end