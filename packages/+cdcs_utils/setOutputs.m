function [x,y,z,opts] = setOutputs(X,Y,Z,others,Kold,cd,Ech,chstuff,opts)

% CDCS/packages/+cdcs_utils/SETOUTPUTS.m
%
% Set output for CDCS. Call the correct routine according to the method 
% specified in opts.solver.

switch lower(opts.solver)
   
    case {'primal', 'dual'}
        [x,y,z,opts] = cdcs_pd.setOutputs(X,Y,Z,others,Kold,cd,Ech,chstuff,opts);
        
    case {'hsde'}
        error('Homogeneous self-dual embedding solver coming soon!')
        
    otherwise
        error('Unknown value for ''options.solver''.')
            
end