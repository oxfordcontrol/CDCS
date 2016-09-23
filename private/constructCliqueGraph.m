%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [adjacencyMatrixC,noOfEdges,edgeCostVectC,incidenceMatrixC] ... 
    = constructCliqueGraph(clique)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparseCoLO 
% Copyright (C) 2009 
% Masakazu Kojima Group
% Department of Mathematical and Computing Sciences
% Tokyo Institute of Technology
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% adjacensy matrix ---> 
% startingTime = cputime; 
% adjacencyMatrixC = sparse(clique.NoC,clique.NoC);
% for p=1:clique.NoC-1
%     for q=p+1:clique.NoC
%         tempSet = intersect(clique.Set{p},clique.Set{q});
%         if ~isempty(tempSet) 
%             adjacencyMatrixC(p,q) = length(tempSet);
%         end
%     end
% end
cliqueMat = sparse(clique.NoC,clique.NoC);
for p=1:clique.NoC
    cliqueMat(p,clique.Set{p}) = 1;
end
% idx = 0; 
% for p=1:clique.NoC
%     sDimE = clique.NoElem(p);
%     cliqueMat(p,clique.Elem(idx+(1:sDimE))) = 1;
%     idx + sDimE; 
% end
adjacencyMatrixC = triu(cliqueMat*cliqueMat',1); 
% fprintf('Part A: cputime = %6.2f\n',cputime-startingTime); 
% <--- adjacensy matrix
% edgeCostVectC --->
% startingTime = cputime; 
edgeCostVectC = reshape(adjacencyMatrixC',1,clique.NoC*clique.NoC); 
nzIdx = find(edgeCostVectC); 
edgeCostVectC = full(edgeCostVectC(nzIdx)); 
% fprintf('Part B: cputime = %6.2f\n',cputime-startingTime); 
%edgeCostVectC
% <--- edgeCostVectC
% the number of edges ---> 
noOfEdges = nnz(adjacencyMatrixC); 
% <--- the number of edges
% incidence matrix ---> 
% startingTime = cputime; 
% incidenceMatrixC = []; 
% edgePointer = 0;
% identityMat = speye(clique.NoC,clique.NoC); 
% for p=1:clique.NoC-1
%     nzIdx = find(adjacencyMatrixC(p,:) > 0);
%     if ~isempty(nzIdx)
%         ll = length(nzIdx); 
%         addMatrix = sparse(clique.NoC,ll);
%         addMatrix(p,:) = 1;        
%         addMatrix = addMatrix - identityMat(:,nzIdx); 
%         incidenceMatrixC = [incidenceMatrixC,addMatrix]; 
%     end
% end
idxNz = find(adjacencyMatrixC'); 
idx1 = ceil(idxNz/clique.NoC)'; 
idx2 = mod(idxNz,clique.NoC)';
idxZero = find(idx2==0);
idx2(idxZero) = clique.NoC;
identityMat = speye(clique.NoC,clique.NoC);
incidenceMatrixC = identityMat(:,idx1) - identityMat(:,idx2); 
% fprintf('Part C: cputime = %6.2f\n',cputime-startingTime); 
% <--- incidentce matrix 
% Ordering columns according to edge costs --->
[edgeCostVectC,permutation] = sort(edgeCostVectC,'descend');
incidenceMatrixC = incidenceMatrixC(:,permutation); 
% 
% edgeCostVectC
% 
% XXXX
% <--- Ordering columns according to edge costs
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
