function [header,myline1,myline2] = printHeader

% CDCS/packages/cdcs_pd/PRINTHEADER
% Print header for iterations

% Set stuff
myline1 = [repmat('=',1,86),'\n'];
myline2 = [repmat('-',1,86),'\n'];
header  = [' iter |   pres   |   dres   |   pcost   |   dcost   |   gap    ',...
    '|   rho    | time (s) |\n'];
