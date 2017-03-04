function [header,myline1,myline2] = printHeader(opts)

% CDCS/packages/+cdcs_utils/PRINTHEADER.m
%
% Print header for solver

switch lower(opts.solver)
   
    case {'primal', 'dual'}
        [header,myline1,myline2] = cdcs_pd.printHeader;
        
    case {'hsde'}
        [header,myline1,myline2] = cdcs_hsde.printHeader;
        
    case {'sos'}
        [header,myline1,myline2] = cdcs_sos.printHeader;
        
    otherwise
        error('Unknown value for ''options.solver''.')
            
end

