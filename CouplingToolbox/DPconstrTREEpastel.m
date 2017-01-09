function pS = DPconstrTREEpastel(pS,Wth,treeMetric)

%This function constructs a kd tree of a pointset of multiple concatenated  
%trials using PASTEL
%(http://kaba.hilvi.org/pastel/pastel/geometrymatlab/pastelgeometrymatlab.htm)


%Inputs:
%pS: (point structure) is a structure with fields:
    %X: a pointset of Nt points, D dimensions and Ntr trials in (points x
    %   dimensions x trials) form
    %Nt: number of points (time indexes) of X
    %D: dimensionality of X
    %Ntr: number of trials of X
%Wth: Theiler window in number of points
%treeMetric: (default) distance metric, either 'euclidean' or 'maximum' norm

%Inputs:
%pS: (point structure) is a structure with fields:
    %X: a pointset of Nt points, D dimensions and Ntr trials in (points x
    %   dimensions x trials) form
    %Nt: number of points (time indexes) of X
    %D: dimensionality of X
    %Ntr: number of trials of X

%Outputs
%pS: (tree structure) a structure with fields
    %X: the pointset of concatenated trials
    %Nt: number of points (time indexes) of the original pointset X
    %D: dimensionality of X
    %Ntr: number of trials of X
    %kdTree: the kd tree structure
    %time2treeIND: a look-up table between time indexes and indexes of pointset of 
    %   concatenated trials 
    %tree2timeIND: a look-up table between indexes of pointset of 
    %   concatenated trials and time indexes  
    %Wth: Theiler window in number of points
    %treeMetric: (default) distance metric, either 'euclidian' or 'maximum' norm 
    

%Permute x so that it is in the (points x trials x dimensions) form
pS.X = permute(pS.X,[1,3,2]);

%Concatenate trials, i.e. reshape X to take the form 
%(points*Ntr x dimensions):
pS.X = reshape(pS.X,[pS.Nt*pS.Ntr,pS.D]);

%Permute x so that it is in the (dimensions x points*trials) form
pS.X = permute(pS.X,[2,1]);


% Construct a kd-tree.
%pS.kdTree = pointkdtree_construct(pS.D);
pS.kdTree = pastelgeometrymatlab('pointkdtree_construct', pS.D);

% Insert the points into the kd-tree.
%idSet = pointkdtree_insert(pS.kdTree, pS.X);
idSet = pastelgeometrymatlab('pointkdtree_insert', pS.kdTree, pS.X);

%Permute x so that it is in the (points*trials x dimensions) form
pS.X = permute(pS.X,[2,1]);

% Refine the subdivision of the kd-tree.
%pointkdtree_refine(pS.kdTree);
pastelgeometrymatlab('pointkdtree_refine',pS.kdTree, 8);

%Construct time-to-tree indexes'look-up table:
IND = 1:pS.Nt;
pS.time2treeIND = nan(pS.Nt,pS.Ntr);
for iTr = 1:pS.Ntr;
    pS.time2treeIND(:,iTr) = IND + (iTr-1)*pS.Nt ;
end


%Construct tree-to-time indexes'look-up table:
pS.tree2timeIND = repmat(IND.',[pS.Ntr,1]);


pS.Wth=Wth;
pS.treeMetric = treeMetric;

