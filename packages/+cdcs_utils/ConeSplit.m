%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clique] = ConeSplit(alpha)
%
% Split a big matrix into multple blocks
% X: any matrix of dimension N #times N
% alpha: a partition -- k1, k2, ..., kn
% the structures of clique are 
%     NoC, maxC, minC, Elem, NoElem and idxMatrix.
%
% Elem is a row vector, which contains all cliques.
% NoElem is NOC-dimensional row vector, which has each size of cliques
%
% For example, we consider three cliques
% c1 ={1,2}, c2 = {1,3,4}, c3 ={5}.
%
% Then 
% clique.Elem =[1,2,1,3,4,5]
% clique.NoElem = [2,3,1]
%


n      = length(alpha);

if n == 1
    clique.NoC = 1;
    clique.Elem = 1:alpha;
    clique.NoElem = alpha;     
    clique.maxC = alpha;
    clique.minC = alpha;
    clique.idxMatrix = ones(alpha);  
else
    Num = nchoosek(n,2);
    Tn  = cumsum([1;alpha(:)]);
    clique.NoC = Num;
    clique.Elem = [];
    clique.NoElem = zeros(Num,1);
    iter = 1;
    for i = 1:n
        for j = i+1:n
            tmp1 = Tn(i):Tn(i+1)-1;
            tmp2 = Tn(j):Tn(j+1)-1;
            clique.Elem = [clique.Elem; tmp1(:);tmp2(:)];
            clique.NoElem(iter) = alpha(i) + alpha(j);
            iter = iter + 1;
        end
    end

    clique.maxC = max(clique.NoElem);
    clique.minC = min(clique.NoElem);
    clique.idxMatrix = ones(sum(alpha));

end