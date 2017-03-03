function [At,C,Ats,Cs] = svecData(At,C,K)

% SVECDATA.m
%
% Transform data from vectorize format to svectorized format for semidefinite
% variables according to the cone K.

%--------------------------------------------
% Some variables
%--------------------------------------------
nSDP = sum(K.s.^2);
nonSDP = K.f + K.l + sum(K.q);

if nSDP > 0
    %--------------------------------------------
    % Build matrix Q such that svec(X)=Q*vec(X).
    % Consider each block separately.
    %--------------------------------------------
    nBlk = length(K.s);
    isr2 = 1./sqrt(2);
    I = zeros(nSDP,1);
    J = zeros(nSDP,1);
    V = zeros(nSDP,1);
    count = 0;
    bdrow = 0;
    bdcol = 0;
    
    % Loop over blocks
    for blk = 1:nBlk
        n = K.s(blk);
        rowcnt = 0;
        
        % Loop over columns
        for j=1:n
            
            % Find indices
            nels = n-j+1;
            rws = repmat(1:nels,2,1);
            rws = rws(:);
            cls = [n*(j-1)+(j+1:n); n*((j+1:n)-1)+j];
            I(count+1:count+2*nels-1) = rowcnt + rws(2:end) + bdrow;
            J(count+1:count+2*nels-1) = [n*(j-1)+j; cls(:)] + bdcol;
            V(count+1:count+2*nels-1) = [1; repmat(isr2,2*nels-2,1)];
            
            % update counters
            rowcnt = rowcnt + nels;
            count = count + 2*nels-1;
        end
        
        % Update block diagonal counter
        bdrow = bdrow + n*(n+1)/2;
        bdcol = bdcol + n^2;
    end
    
    % Build sparse Q
    Q  = sparse(I,J,V,bdrow,nSDP);
    
    % svec the data
    Ats = Q*At(nonSDP+1:end,:);
    Cs  = Q*C(nonSDP+1:end,:);
    
else
    Ats = [];
    Cs = [];
    
end

%--------------------------------------------
% Modify At and c
%--------------------------------------------
At  = [At(1:nonSDP,:); Ats];
C   = [ C(1:nonSDP); Cs];


end