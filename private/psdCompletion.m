function x = psdCompletion(x,K,clique)

% PSDCOMPLETION.m
%
% This is a modified version of SparseCoLO/psdCompletion, adapted for use in
% CDCS.

debugSW = 0;            % if 1, print diagnostics
epsilon = 1.0e-6;       % perturbation to ensure that XMat(U,U) is positive definite
% If
% ??? Error using ==> chol
% Matrix must be positive definite.
%
% Error in ==> psdCompletion at 61
%                         LMat = chol(XMat(U,U)+epsilon*speye(nDim,nDim));
%
% Then take a larger epsilon!


rowPointer = K.f+K.l+sum(K.q);
uVect = x(rowPointer+1:end,:);      % the SDP portion
x = x(1:rowPointer);                % the non-SDP part
if isfield(K,'s') && ~isempty(K.s) && sum(K.s)>0
    noOfSDPcones = length(K.s);
    rowPointer = 0;
    for kk=1:noOfSDPcones
        sDim = K.s(kk);
        % Code from <New version> of SparseCoLO/psdCompletion.m
        if clique{kk}.NoC > 1
            XMat = full(reshape(uVect(rowPointer+1:rowPointer+sDim*sDim,1),sDim,sDim));
            if ~isfield(clique{kk},'NoCliqueInForest')
                %                1
                [adjacencyMatrixC,noOfEdges,edgeCostVectC,incidenceMatrixC] ...
                    = constructCliqueGraph(clique{kk});
                %
                % We assume that the clique graph is connected
                %
                randSeed = 2009;
                [treeValue,adjacencyMatrixT,edgeCostVectT,incidenceMatrixT] ...
                    = maxSpanningTree(clique{kk},adjacencyMatrixC,edgeCostVectC,incidenceMatrixC,randSeed);
                kDim = size(incidenceMatrixT,2);
                %
                %
                %
                for i=1:clique{kk}.NoC-1
                    edgeSet = find(incidenceMatrixT(i,:) ~= 0);
                    if ~isempty(edgeSet)
                        j = edgeSet(1);
                        I = find(incidenceMatrixT(:,j)' ~= 0);
                        i2 = I(2);
                        U = intersect(clique{kk}.Set{i},clique{kk}.Set{i2});
                        S = setdiff(clique{kk}.Set{i},U);
                        T = setdiff(clique{kk}.Set{i2},U);
                        nDim = length(U);
                        % epsilon = 1.0e-10;
                        LMat = chol(XMat(U,U)+epsilon*speye(nDim,nDim));
                        XMat(S,T) = (XMat(S,U)/LMat)*(LMat'\XMat(U,T));
                        XMat(T,S) = XMat(S,T)';
                        clique{kk}.Set{i2} = union(clique{kk}.Set{i},clique{kk}.Set{i2});
                        incidenceMatrixT(i2,:) = incidenceMatrixT(i2,:) + incidenceMatrixT(i,:);
                        incidenceMatrixT(i,:) = sparse(1,kDim);
                    end
                end
            else
                %                2
                fPointer = 0;
                for ii=1:length(clique{kk}.NoCliqueInForest)
                    if clique{kk}.NoCliqueInForest(ii) > 1
                        tempClique.NoC = 0;
                        tempClique.NoElem = [];
                        tempClique.maxC = 0;
                        tempClique.minC = 1.0e10;
                        setIdx = 0;
                        for j=fPointer+1:fPointer+clique{kk}.NoCliqueInForest(ii)
                            tempClique.NoC = tempClique.NoC + 1;
                            tempClique.NoElem = [tempClique.NoElem,length(clique{kk}.Set{j})];
                            setIdx = setIdx + 1;
                            tempClique.Set{setIdx} = clique{kk}.Set{j};
                        end
                        tempClique.maxC = max(tempClique.NoElem);
                        tempClique.minC = min(tempClique.NoElem);
                        [adjacencyMatrixC,noOfEdges,edgeCostVectC,incidenceMatrixC] ...
                            = consructCliqueGraph(tempClique);
                        randSeed = 2009;
                        [treeValue,adjacencyMatrixT,edgeCostVectT,incidenceMatrixT] ...
                            = maxSpanningTree(tempClique,adjacencyMatrixC,edgeCostVectC,incidenceMatrixC,randSeed);
                        kDim = size(incidenceMatrixT,2);
                        for i=1:tempClique.NoC-1
                            edgeSet = find(incidenceMatrixT(i,:) ~= 0);
                            if ~isempty(edgeSet)
                                j = edgeSet(1);
                                I = find(incidenceMatrixT(:,j)' ~= 0);
                                i2 = I(2);
                                U = intersect(tempClique.Set{i},tempClique.Set{i2});
                                S = setdiff(tempClique.Set{i},U);
                                T = setdiff(tempClique.Set{i2},U);
                                nDim = length(U);
                                LMat = chol(XMat(U,U)+epsilon*speye(nDim,nDim));
                                XMat(S,T) = (XMat(S,U)/LMat)*(LMat'\XMat(U,T));
                                XMat(T,S) = XMat(S,T)';
                                tempClique.Set{i2} = union(tempClique.Set{i},tempClique.Set{i2});
                                incidenceMatrixT(i2,:) = incidenceMatrixT(i2,:) + incidenceMatrixT(i,:);
                                incidenceMatrixT(i,:) = sparse(1,kDim);
                            end
                        end
                        clear tempClique
                    end
                    fPointer = fPointer + clique{kk}.NoCliqueInForest(ii);
                end
                
            end
            
            % Debug?
            if debugSW == 1
                d = eig(XMat);
                fprintf('the minimum eigenvalue of a completed SDP variable matrix = %+6.1e\n',full(min(d')));
            end
            
            % Set
            x = [x; reshape(XMat,sDim*sDim,1)];
        else
            x = [x; uVect(rowPointer+1:rowPointer+sDim*sDim,1)];
        end
        rowPointer = rowPointer + sDim*sDim;
    end
end

% END FUNCTION
end
