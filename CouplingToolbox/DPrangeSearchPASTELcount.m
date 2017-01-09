function [nncount] = DPrangeSearchPASTELcount(pS, r, L, queryIND, treeMetric)

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


if nargin<5
    treeMetric = pS.treeMetric;
end

if strcmpi(treeMetric,'euclidean')
    r = r.^2;
end

%Preallocate memory
nncount = zeros(L,pS.Ntr);

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
        
        %...search the number of neighbors within distance r...
        %  countSet = count_nearest(self, querySet, maxDistanceSet, varargin)
        %varargin = {'norm'}
        % Optional input arguments in 'key'-value pairs
        %nncount(iT,iTr) = range_search(pS.kdTree, pS.time2treeIND(queryIND(iT),iTr), r,treeMetric);
        nncount(iT,iTr) = pastelgeometrymatlab('pointkdtree_count_nearest', pS.kdTree, pS.time2treeIND(queryIND(iT),iTr), r, treeMetric);
        
        %...show back the Theiler zone points...
        %pointkdtree_show(kdTree, idSet)
        %pointkdtree_show(pS.kdTree, pS.time2treeIND(theilerZone,iTr))
        pastelgeometrymatlab('pointkdtree_show', pS.kdTree, pS.time2treeIND(theilerZone,iTr));
   
    end
end


