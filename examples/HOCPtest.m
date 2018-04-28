
% HOCP TEST
%
% Run some example with full pattern

% Preliminaries
clc;
opts.maxIter = 1e+3;
opts.relTol  = 1e-3;





% ---------------------------------------------------------------------------- %
%                       Conic program with banded SDP
% ---------------------------------------------------------------------------- %
% Parameters
m   = 30;                      % # constraints
K.s = [1000];                    % PSD cones
bandWidth = [K.s-1];              % bandWidth for SDP cones


% Setup
fprintf('\nSetting up random conic problem with banded SDP cones, m=%i...',m);
tsetup = tic;
[At,b,c,K] = bandedSDP(m,K,bandWidth);
tsetup = toc(tsetup);
fprintf('done in %.2f seconds. \n',tsetup);

% solution by admm
opts.solver = 'primal';
[x1,y1,z1,info1] = cdcs(At,b,c,K,opts);

opts.solver = 'hocp';      % using higher oder conic programs
opts.NoP = 10;             % 10 blocks
[x2,y2,z2,info2] = cdcs(At,b,c,K,opts);

%  The following setting is using second order cone programs  -- not
%  efficient
%  opts.NoP = K.s;            % second order cone programs
%  [x2,y2,z2,info2] = cdcs(At,b,c,K,opts);
