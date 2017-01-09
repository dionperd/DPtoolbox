function [indNNtree distNN nnX indNNtime] = DPnnSearchTSTOOL(pS, k, L, queryIND, treeMetric, maxDistanceSet)

%This function performs fixed mass nearest neighbors search to a kdTree 
%of a pointset of concatenated trials using TSTOOL toolbox
%(http://www.physik3.gwdg.de/tstool/index.html)



%Inputs:
%pS: (tree structure) a structure with fields
    %X: the pointset of concatenated trials ((points+Wth)*trials x
%       dimensions) form
    %Nt: number of points (time indexes) of the original pointset X
    %D: dimensionality of X
    %Ntr: number of trials of X
    %kdTree: the kd tree structure
    %time2treeIND: a look-up table between time indexes and indexes of pointset of 
    %   concatenated trials 
    %tree2timeIND: a look-up table between indexes of pointset of 
    %   concatenated trials and time indexes  
    %Wth: Theiler window in number of points
    %treeMetric: distance metric, either 'euclidian' or 'maximum' norm
    
%k: number of neighbors for the fixed mass
%L: number of query points in time
%queryIND: indexes of query points in time
%treeMetric: distance metric, either 'euclidean' or 'maximum' norm
%maxDistanceSet: a matrix of size (1,L) of maximum range for the search

%Outputs
%indNNtree: the indexes of the nearest neighbors in the tree, a matrix of
%(L,Ntr,k) size
%distNN: the distances of the nearest neighbors, a matrix of (L,k,Ntr) size
%nnX = pointset of nearest neighbors, a matrix of (L,D,Ntr,k) size
%indNNtime: the indexes of the nearest neighbors in time, a matrix of (L,Ntr,k)
%size


%Preallocate memory
indNNtree = zeros(pS.Ntr,k,L);
if nargout>1
    distNN = zeros(pS.Ntr,k,L);
    if nargout>2
        nnX = nan(k,pS.D,pS.Ntr,L);
        if nargout>3
            indNNtime = zeros(pS.Ntr,k,L);
        end
    end
end


%For each query time point...
for iT = 1:L;    
    
     %...search the k nearest neighbors...
     % [index, distance] = nn_search(pointset, atria, query_indices, k, exclude, epsilon);
     [indNNtree(:,:,iT), distNN(:,:,iT)] = nn_search(pS.X, pS.kdTree, pS.time2treeIND(queryIND(iT),:).', k, pS.Wth, 0);
     
     if nargout>2
         %...and for each trial...
         for iTr = 1:pS.Ntr;
             %...take only the non zero indices...
             tempNOzeroINDs = indNNtree(iTr,:,iT)>0;
             %...get the coordinates of the  k nearest neighbors...
             nnX(tempNOzeroINDs,:,iTr,iT) = pS.X(indNNtree(iTr,tempNOzeroINDs,iT),:);
             if nargout>3
                 %...convert the tree indices of the k nearest neighbors into time indices...
                 indNNtime(iTr,tempNOzeroINDs,iT) = pS.tree2timeIND(indNNtree(iTr,tempNOzeroINDs,iT));
             end
             
         end
     end
     
end

%Permute to take the desired form
indNNtree=permute(indNNtree,[3,1,2]);
if nargout>1
    distNN=permute(distNN,[3,1,2]);
    if nargout>2
        nnX = permute(nnX,[4,2,3,1]);
        if nargout>3
            indNNtime=permute(indNNtime,[3,1,2]);
        end
    end
end

