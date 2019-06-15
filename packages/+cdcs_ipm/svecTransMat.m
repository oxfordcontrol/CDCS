function Q = svecTransMat(n)

%% SVECTRANSMAT.m Matrix to turn vec(X) into svec(X)
%
% Compute transformation matrix Q so that for a given n\times n symmetric matrix
% X one has
% 
%   svec(X) = Q*vec(X)

isr2 = 1./sqrt(2);
I = zeros(n^2,1); 
J = zeros(n^2,1);
V = zeros(n^2,1);
count = 0;
rows = 0;
for j=1:n                   % for each column of X 
    
    nels = n-j+1;
    r = repmat(1:nels,2,1); r = r(:);
    c = [n*(j-1)+(j+1:n); n*((j+1:n)-1)+j];
    I(count+1:count+2*nels-1) = rows + r(2:end);
    J(count+1:count+2*nels-1) = [n*(j-1)+j; c(:)];
    V(count+1:count+2*nels-1) = [1; repmat(isr2,2*nels-2,1)];
    rows = rows + nels;
    count = count + 2*nels-1;

end
Q = sparse(I,J,V,n*(n+1)/2,n^2);