%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clique] = cliquesFromSpMatD(sparsityPatternMat)

% Modified by M. Kojima,March 25, 2010

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

%
% 2008-06-13 Waki
% Caution!
% I have changed the structures of clique.
% 
% Before:
% the structures of clique are 
%     NoC, maxC, minC, Set and idxMatrix.
%
% After:
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
% If you want to access the i-th clqiue, execute the following command:
%
% idx = sum(clique.NoElem(1:i-1));
% clique.Elem(idx+(1:clique.NoElem(i)));
%
% In this example, if you want to access the second clique, you should execute
%
% idx = clique.NoElem(1);
% clique.Elem(idx+(1:clique.NoElem(2)));
%

% tic

%spones(sparsityPatternMat);
nDim = size(sparsityPatternMat,1);
sparsityPatternMat = spones(sparsityPatternMat) + (2*nDim+1)*speye(nDim);
%%%%%
orderingSW = 0;
if orderingSW == 0%
% pars.edgeRemoveSW == 0;
    %% minimum degree ordering
%    orderingSW = 0;
    I = symamd(sparsityPatternMat);
elseif orderingSW == 1
    %% sparse reverse Cuthill-McKee ordering
%    orderingSW = 1;
    I = symrcm(sparsityPatternMat);   
elseif orderingSW == 3
% pars.edgeRemoveSW = 3;
    %% sparse reverse Cuthill-McKee ordering only for the sensor network
    %% problem
%    orderingSW = 1;
    pars.sDim = 3; 
    I = symrcm(sparsityPatternMat(1:nDim-pars.sDim,1:nDim-pars.sDim));
    I = [I,(nDim-pars.sDim+1):nDim];     
%    I = symamd(sparsityPatternMat);
% if orderingSW == 0
%     
%     %% minimum degree ordering
%     I = symamd(sparsityPatternMat);
% elseif orderingSW == 1
%     %% sparse reverse Cuthill-McKee ordering
%     I = symrcm(sparsityPatternMat);
end
%% cholesky decomposition
[R,p] = chol(sparsityPatternMat(I,I));
if (p > 0)
    error('Correlative sparsity matrix is not positive definite.');
end

% debug = 1;
% if debug == 1
%     figure(1);
%     spy(R); 
% end

%%
%% Step3
%% Finding the maxmal clieques
%%
%% put 1 for nonzero element of R
% if ORIGINAL == 1
%     Cliques = spones(R);
%     [value,orig_idx] = sort(I);
%     remainIdx = 1;
%     for i=2:nDim
%         checkSet = Cliques(i,i:nDim);
%         one = find(checkSet);
%         noOfone = length(one);
%         cliqueResult = Cliques(1:i-1,i:nDim)*checkSet';
%         yesno = find(cliqueResult == noOfone);
%         %%
%         %% Remove the set included into other set.
%         %%
%         if ~any(yesno)
%             remainIdx = [remainIdx;i];
%         end
%     end
%     clique.Set = Cliques(remainIdx,orig_idx);
%     %%
%     %% Clique Information
%     %%
%     clique.NoC  = size(clique.Set,1);
%     sumClique = full(sum(clique.Set,2));
%     clique.maxC = full(max(sumClique));
%     clique.minC = full(min(sumClique));
%     %% Information on indexing variables
%     % Cliques = spones(R+R');
%     % Cliques = Cliques(orig_idx,orig_idx);
%     Cliques = sparse(nDim,nDim);
%     for i=1:clique.NoC
%         idx = find(clique.Set(i,:));
%         sDimE = length(idx);
%         Cliques(idx,idx) = ones(sDimE,sDimE);
%     end
%     Cliques = triu(Cliques);
%     % % spy(Cliques);
%     pointer = 0;
%     for i=1:nDim
%         nnzRowIdx = find(Cliques(i,:));
%         noNnz = length(nnzRowIdx);
%         Cliques(i,nnzRowIdx) = [pointer+1:pointer+noNnz];
%         pointer = pointer + noNnz;
%     end
%     % full(Cliques)
%     clique.idxMatrix = Cliques;
% elseif ORIGINAL == 0
    Cliques = spones(R);
%     
%     full(Cliques)
%     
    [value,orig_idx] = sort(I);
    remainIdx = 1;
    
% fprintf('\n0: %10.5fseconds\n\n',toc);
% tic
    for i=2:nDim
        idx = i:nDim;
        one = find(Cliques(i,idx));
        noOfone = length(one);
        %
        % 2008-06-08 Waki
        % replace multiplication by sum.
        %
        cliqueResult = sum(Cliques(remainIdx,idx(one)),2);
%         
%         full(Cliques(remainIdx,idx(one)))
%         idx
%         full(one)
%         idx(one)
%         noOfone
%         full(cliqueResult)
%         
%         XXXXX
%         
%        if isempty(find(full(cliqueResult) == noOfone,1))
        if isempty(find(cliqueResult == noOfone,1))
            remainIdx = [remainIdx;i];
        end
    end
    
% fprintf('\n1: nDim = %d, %10.5fseconds\n\n',nDim,toc);
% tic
    
%     debugSW = 1;
%     if debugSW == 1
%         figure(2);
%         spy(R(remainIdx,:)); 
%         figure(3);
%         spy(Cliques(remainIdx,:));
%         RR = Cliques(remainIdx,orig_idx)'*Cliques(remainIdx,orig_idx); 
%         figure(4); 
%         spy(triu(RR));
%     end

    cSet = Cliques(remainIdx,orig_idx);
%       cSet = Cliques(remainIdx,:);
%       full(cSet)
%    XXXXX
%       
    %%
    %% Clique Information
    %%
    clique.NoC  = length(remainIdx);
    [I,J] = find(cSet');
    clique.Elem = I;
    clique.NoElem = full(sum(cSet,2));
    clique.NoElem = clique.NoElem';
    clique.maxC = full(max(clique.NoElem));
    clique.minC = full(min(clique.NoElem));

% Modified by M. Kojima,March 25, 2010 ---> 
%     tic; 
%     modifySW = 1;
%     if modifySW == 0
%         Cliques = sparse(nDim,nDim);
%         idx = 0;
%         for i=1:clique.NoC
%             s = clique.NoElem(i);
%             tmp = clique.Elem(idx+(1:s));
%             idx = idx + s;
%             %
%             % 2008-06-08 Waki
%             % replace ones by 1
%             %
%             Cliques(tmp,tmp) = 1;
%         end
%         Cliques = tril(Cliques);
%         
%         %         debugSW = 1;
%         %         if debugSW == 1
%         %             figure(10);
%         %             spy(Cliques);
%         %         end
%     elseif modifySW == 1
        Cliques = tril(Cliques(remainIdx,orig_idx)'*Cliques(remainIdx,orig_idx));
        %         debugSW = 1;
        %         if debugSW == 1
        %             figure(20);
        %             spy(Cliques);
        %         end
%    end
%    fprintf('\n1: modifySW = %d, %10.5f seconds\n\n',modifySW, toc);
% <--- Modified by M. Kojima,March 25, 2010
    
    %
    % 2008-06-10 Waki
    % Replace the part of renumbering Cliques by a faster version.
    %
    [I,J,V] = find(Cliques);
    s = length(V);
    clique.idxMatrix = sparse(J,I,(1:s),nDim,nDim,s);
    
% fprintf('\n2: %10.5fseconds\n\n',toc);
% tic
        
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% 
% if MEMORY == 1
%     d = whos('*');
%     for i=1:length(d)
%         mem = mem - d(i).bytes;
%     end
%     fprintf('$$ The total Memory in cliquesFromSpMatD = %4d KB\n',-mem/1024);
% end

% Listing all cliques by Kojima 2008/07/16 ---> 
clique.Set{1} = clique.Elem(1:clique.NoElem(1))';
for p=2:clique.NoC
    idx = sum(clique.NoElem(1:p-1));
	clique.Set{p} = clique.Elem(idx+(1:clique.NoElem(p)))';
end 
% <--- Listing all cliques by Kojima 2008/07/16
% fprintf('\n3: clique.NoC = %d, %10.5fseconds\n\n',clique.NoC,toc);
% XXXXX

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
