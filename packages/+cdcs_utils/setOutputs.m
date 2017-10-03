function [x,y,z,info,opts] = setOutputs(X,Y,Z,others,Kold,c,Ech,chstuff,info,opts)

% CDCS/packages/+cdcs_utils/SETOUTPUTS.m
%
% Set output for CDCS. Call the correct routine according to the method 
% specified in opts.solver.

switch lower(opts.solver)
   
    case {'primal', 'dual'}
        [x,y,z,info,opts] = cdcs_pd.setOutputs(X,Y,Z,others,Kold,c,Ech,chstuff,info,opts);
        
    case {'hsde'}
        [x,y,z,info,opts] = cdcs_hsde.setOutputs(X,Y,Z,others,Kold,c,Ech,chstuff,info,opts);
        
    case {'sos'}
        [x,y,z,info,opts] = cdcs_sos.setOutputs(X,Y,Z,Kold,info,opts);
        
    otherwise
        error('Unknown value for ''options.solver''.')
            
end