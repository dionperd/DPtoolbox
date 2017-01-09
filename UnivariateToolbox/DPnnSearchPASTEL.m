function [indNNtree distNN nnX indNNtime] = DPnnSearchPASTEL(pS, k, L, queryIND, treeMetric, maxDistanceSet)

%This function performs fixed mass nearest neighbors search to a kdTree 
%of a pointset of concatenated trials using PASTEL
%(http://kaba.hilvi.org/pastel/pastel/geometrymatlab/pastelgeometrymatlab.htm)



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
    %treeMetric: (default) distance metric, either 'euclidean' or 'maximum'
    %norm 
%k: number of neighbors for the fixed mass
%L: number of query points in time
%queryIND: indexes of query points in time
%treeMetric: distance metric, either 'euclidean' or 'maximum' norm
%maxDistanceSet: a matrix of size (1,L) of maximum range for the search

%Outputs
%indNNtree: the indexes of the nearest neighbors in the tree, a matrix of (L,Ntr,k) size
%indNNtime: the indexes of the nearest neighbors in time, a matrix of (L,Ntr,k)
%size
%distNN: the distances of the nearest neighbors, a matrix of (L,k,Ntr) size
%nnX = pointset of nearest neighbors, a matrix of (L,D,Ntr,k) size



if nargin<6
    maxDistanceSet = Inf;
end

if nargin<5
    treeMetric = pS.treeMetric;
end
    
if (maxDistanceSet ~= Inf) && strcmpi(treeMetric,'euclidean')
    maxDistanceSet = maxDistanceSet^2;
end

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
    
     %...construct the Theiler zone in time...
     theilerZone = queryIND(iT)-pS.Wth:queryIND(iT)+pS.Wth;
     %...by getting rid of no positive and great than Nt indexes...
     theilerZone=theilerZone( (theilerZone>0) & (theilerZone<=pS.Nt));
     
     
     %...and for each trial...
     for iTr = 1:pS.Ntr;
         
         %...hide the Theiler zone points...
         %pointkdtree_hide(kdTree, idSet)
         %pointkdtree_hide(kdTree, pS.time2treeIND(theilerZone,iTr))
         pastelgeometrymatlab('pointkdtree_hide', pS.kdTree, int64(pS.time2treeIND(theilerZone,iTr)));
         
         %...search the k nearest neighbors...
         %[neighborSet, distanceSet] = search_nearest(self, querySet, varargin)
         %varargin = {kNearest, maxDistanceSet, norm}
         % Optional input arguments in 'key'-value pairs
         %[indNNtree(iTr,:,iT) distNN(iTr,:,iT)]  = pointkdtree_search_nearest(pS.kdTree, pS.time2treeIND(queryIND(iT),iTr), 'kNearest',k, 'maxDistanceSet', maxDistanceSet,'norm',treeMetric);
         %[indNNtree(iTr,:,iT) distNN(iTr,:,iT)]  = pointkdtree_search_nearest(pS.kdTree, pS.time2treeIND(queryIND(iT),iTr),  maxDistanceSet,k,treeMetric);
         [indNNtree(iTr,:,iT) distNN(iTr,:,iT)] = pastelgeometrymatlab('pointkdtree_search_nearest', pS.kdTree, int64(pS.time2treeIND(queryIND(iT),iTr)), maxDistanceSet, int64(k), treeMetric);
         
         %...show back the Theiler zone points...
         %pointkdtree_show(kdTree, idSet)
         %pointkdtree_show(pS.kdTree, pS.time2treeIND(theilerZone,iTr))
         pastelgeometrymatlab('pointkdtree_show', pS.kdTree, int64(pS.time2treeIND(theilerZone,iTr)));
         
         %...take only the non zero indices...
         tempNOzeroINDs = indNNtree(iTr,:,iT)>0;
         
         if nargout>2
             
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




