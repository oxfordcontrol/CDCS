function [header,myline1,myline2] = printHeader(opts)

% CDCS/packages/+cdcs_utils/PRINTHEADER.m
%
% Print header for solver

switch lower(opts.solver)
   
    case {'primal', 'dual'}
        [header,myline1,myline2] = cdcs_pd.printHeader;
        
    case {'hsde'}
        %error('Homogeneous self-dual embedding solver coming soon!')
        [header,myline1,myline2] = cdcs_hsde.printHeader;
        
    otherwise
        error('Unknown value for ''options.solver''.')
            
end

