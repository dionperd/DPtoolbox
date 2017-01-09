function [nncount indNNtree distNN nnX indNNtime] = DPrangeSearchPASTEL(pS, r, L, queryIND, treeMetric)

%This function performs fixed range nearest neighbors search to a kdTree 
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
    %treeMetric: distance metric, either 'euclidean' or 'maximum' norm   
%r: range of search
%L: number of query points in time
%queryIND: indexes of query points in time
%treeMetric: (default) distance metric, either 'euclidean' or 'maximum' norm


%Outputs
%if k is the number of points for each search per query point,
%then
%nncount:  a (L,Ntr) matrix of k, per search
%indNNtree: the indexes of the within range neighbors in the tree, (L,Ntr)
%   cell of matrices of (1,k) size
%distNN: the distances of the within range neighbors, a matrix of (L,k,Ntr) size
%nnX = pointsets of nearest neighbors, (L,Ntr) cell of matrices of (k,D)
%   size
%indNNtime: the indexes of the within range neighbors in time, (L,Ntr) cell
%   of matrices of (1,k) size


if nargin<5
    treeMetric = pS.treeMetric;
end

if strcmpi(treeMetric,'euclidean')
    r = r.^2;
end

%Preallocate memory
nncount = zeros(L,pS.Ntr);
if nargout>1
    indNNtree = cell(L,pS.Ntr);
    if nargout>2
        distNN = cell(L,pS.Ntr);
        if nargout>3
            nnX = cell(L,pS.Ntr);
            if nargout>4
                indNNtime = cell(L,pS.Ntr);
            end
        end
    end
end

%For each query time point...
for iT = 1:L;
    
    %...construct the Theiler zone in time...
    theilerZone = queryIND(iT)-pS.Wth:queryIND(iT)+pS.Wth;
    %...by getting rid of no positive and great than Nt indexes...
    theilerZone=theilerZone( (theilerZone>0) & (theilerZone<=pS.Nt));
    
    for iTr = 1:pS.Ntr;
        
        %...hide the Theiler zone points...
        %pointkdtree_hide(kdTree, idSet)
        %pointkdtree_hide(kdTree, pS.time2treeIND(theilerZone,iTr))
        pastelgeometrymatlab('pointkdtree_hide', pS.kdTree, int64(pS.time2treeIND(theilerZone,iTr)));
        
        %...search the k nearest neighbors within distance r...
        [indNNtree{iT,iTr} tempDist] = pastelgeometrymatlab('pointkdtree_search_nearest', pS.kdTree, int64(pS.time2treeIND(queryIND(iT),iTr)), r, int64(L*pS.Ntr), treeMetric);
        
        %...get only the neigbors with non zero indices...
        noZeroINDs = indNNtree{iT,iTr}>0;
        %...return the number of within range neihbors...
        nncount(iT,iTr) = sum(noZeroINDs);
        
        %...show back the Theiler zone points...
        %pointkdtree_show(kdTree, idSet)
        %pointkdtree_show(pS.kdTree, pS.time2treeIND(theilerZone,iTr))
        pastelgeometrymatlab('pointkdtree_show', pS.kdTree, int64(pS.time2treeIND(theilerZone,iTr)));
        
        if nargout>1
            indNNtree{iT,iTr} = indNNtree{iT,iTr}(noZeroINDs);
            
            if nargout>2
                distNN{iT,iTr}=tempDist(noZeroINDs);
                
                if nargout>3
                    %...get the coordinates of the  k nearest neighbors...
                    nnX{iT,iTr} = pS.X(indNNtree{iT,iTr},:);
                    
                    if nargout>4
                        %...convert the tree indices of the k nearest neighbors into time indices...
                        indNNtime{iT,iTr}  = pS.tree2timeIND(indNNtree{iT,iTr});
                    end
                end
                
            end
        end
    end
    
end
end


