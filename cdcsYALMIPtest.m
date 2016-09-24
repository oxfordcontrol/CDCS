function cdcsYALMIPtest

% CDCSTALMIPTEST
%
% Test calling CDCS from YALMIP
% Solve a simple QP
%
% min (1/2)x'Q x 
% subject to Ax = b
%
% with Q symmetric with Q>=0. The optimal solution x* can obviously be found 
% with the KKT system
%
% [Q A'] [x] = [0]
% [A 0 ] [y]   [b]
%
% so we check that CDCS runs succesfully by comparing the solution to x*.

% Setup
clear; clc; yalmip ('clear');
opts = sdpsettings('solver','CDCS','verbose',1);

% Create the problem data
n = 100;
m = max(floor(n/3),1);
A = rand(m,n);
b = rand(m,1);
Q = rand(n); Q = Q+Q'; Q = Q - (min(eig(Q)) - 1)*eye(n);

% Find the analytical solution
M = Q\(A.');
y = (A*M)\b;
xa = M*y;

% Setup and solve with YALMIP using CDCS
x = sdpvar(n,1);
optimize(A*x==b,x'*Q*x/2,opts);
xv = value(x);

% Display solution
fprintf('Analytical objective  : %.4e\n',xa'*Q*xa/2);
fprintf('CDCS objective        : %.4e\n',xv'*Q*xv/2);
fprintf('||x* - x_cdcs ||_inf: %.4e\n',norm(xa-xv,'inf'))

end
