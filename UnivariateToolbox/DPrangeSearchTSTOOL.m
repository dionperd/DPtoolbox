function [nncount indNNtree distNN nnX indNNtime] = DPrangeSearchTSTOOL(pS, r, L, queryIND, treeMetric)

%This function performs fixed range nearest neighbors search to a kdTree 
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
    %treeMetric: distance metric, either 'euclidean' or 'maximum' norm   
%r: range of search
%L: number of query points in time
%queryIND: indexes of query points in time
%treeMetric: distance metric, either 'euclidian' or 'maximum' norm


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
    
    %...search the number of neighbors within distance r...
    %  [count, neighbors] = range_search(pointset, atria, query_indices, r, exclude)
    [nncount(iT,:), neighbors] = range_search(pS.X, pS.kdTree, pS.time2treeIND(queryIND(iT),:).', r, pS.Wth);
    
    if nargout>1
        
        %...and for each trial...
        for iTr = 1:pS.Ntr;
            
            %...get the distances of the within range neighbors after you sort them in increasing order...
            [tempDist sortedIND]= sort(neighbors{iTr,2});
            
            %...get the tree/pointset (sorted) indices of the within range neighbors...
            indNNtree{iT,iTr} = neighbors{iTr,1}(sortedIND);
            
            if nargout>2
                distNN{iT,iTr}=tempDist;
                
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


