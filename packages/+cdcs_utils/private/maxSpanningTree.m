%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [treeValue,adjacencyMatrixT,edgeCostVectT,incidenceMatrixT,basisIdx,BInv] ... 
    = maxSpanningTree(clique,adjacencyMatrixC,edgeCostVectC,incidenceMatrixC,randSeed)

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

rand('state',randSeed);
edgeCostVect = edgeCostVectC;
noOfEdges = length(edgeCostVect); 
BInv = speye(clique.NoC-1,clique.NoC-1); 
edgeDegree = zeros(1,noOfEdges); 
edgeDegreeMinus = zeros(1,noOfEdges); 
nodeDegree = zeros(clique.NoC,1); 
noOfTreeEdges = 0; 
basisIdx = zeros(clique.NoC-1,1);
dterminedRows = []; % find(basisIdx'); 
treeValue = 0;
controlSW = 1;
edgeCostPointer = 0; 
iteration = 0;

while controlSW == 1
    iteration = iteration + 1;
    if isempty(edgeCostVect)
%         fprintf('## iteration = %d at maxSpanningTree4\n',iteration); 
        error('## the conversion fails because the clique graph is not connected!');
    end
    currentEdgeCost = edgeCostVect(1);
    j = find(edgeCostVect < currentEdgeCost,1);
    if isempty(j)
        maxIdx = [edgeCostPointer+1:noOfEdges];
        edgeCostPointer = [];
        edgeCostVect = [];
    else
        maxIdx = [edgeCostPointer+1:edgeCostPointer+j-1];
        edgeCostVect = edgeCostVect(j:noOfEdges-edgeCostPointer);
        edgeCostPointer = edgeCostPointer+j-1;
    end
    while (~isempty(maxIdx)) && (controlSW == 1) 
        ll = length(maxIdx); 
        if ll == 1
            i = maxIdx(1); 
            maxIdx = [];            
        else
            [value,j] = min(edgeDegree(maxIdx)); 
            i = maxIdx(j); 
            maxIdx = setdiff(maxIdx,[i]); 
        end
        pivotCol = BInv * incidenceMatrixC(1:clique.NoC-1,i); 
        pivRowIdx = find(pivotCol' ~= 0); 
        pivRowIdx = setdiff(pivRowIdx,dterminedRows); 
        if ~isempty(pivRowIdx)
            ll = length(pivRowIdx);
            if ll == 1
                q = pivRowIdx(1);
            else
                [value,j] = max(rand(1,ll));
                q = pivRowIdx(j);
            end
            basisIdx(q,1) = i;
            dterminedRows = find(basisIdx');
            addDegreeIdx = find(incidenceMatrixC(:,i)');
            nodeDegree(addDegreeIdx) = nodeDegree(addDegreeIdx) + 1;
            idxPlus = find(incidenceMatrixC(:,i)' == 1);
            idxMinus = find(incidenceMatrixC(:,i)' == -1);
            edgeDegree = edgeDegree + abs(incidenceMatrixC(idxPlus,:)); 
            edgeDegree = edgeDegree + abs(incidenceMatrixC(idxMinus,:)); 
            treeValue = treeValue + edgeCostVectC(i);
            pivotMatrix = speye(clique.NoC-1,clique.NoC-1);
            pivotMatrix(:,q) = -pivotCol/pivotCol(q,1);
            pivotMatrix(q,q) = 1/pivotCol(q,1);
            BInv = pivotMatrix*BInv; 
        end
        if length(dterminedRows) == clique.NoC-1
            controlSW = 0;
        end
    end
end

edgeCostVectT = edgeCostVectC(1,basisIdx'); 
incidenceMatrixT = incidenceMatrixC(:,basisIdx'); 

adjacencyMatrixT=sparse(clique.NoC,clique.NoC); 
for p=1:clique.NoC-1
    idx = find(incidenceMatrixT(:,p)'); 
    adjacencyMatrixT(idx(1),idx(2)) = edgeCostVectT(p);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

