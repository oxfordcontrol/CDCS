function [At,b,c,K,Ech,chstuff,opts] = preprocess(At,b,c,K,opts)

% CDCS/packages/+cdcs_utils/PREPROCESS.m
%
% Preprocess data: chordalize and rescale. Call functions according to their
% method (homogeneous self-dual method uses different rescaling!)

switch lower(opts.solver)
   
    case {'primal', 'dual'}
        [At,b,c,K,Ech,chstuff,opts] = cdcs_pd.preprocess(At,b,c,K,opts);
        
    case {'hsde'}
        [At,b,c,K,Ech,chstuff,opts] = cdcs_hsde.preprocess(At,b,c,K,opts);
        
    case {'sos'}
        [At,b,c,K,Ech,chstuff,opts] = cdcs_sos.preprocess(At,b,c,K,opts);
        
    otherwise
        error('Unknown value for ''options.solver''.')           
end

