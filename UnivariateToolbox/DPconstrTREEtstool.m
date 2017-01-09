function pS = DPconstrTREEtstool(pS,Wth,treeMetric)

%This function constructs a kd tree of a pointset of multiple concatenated  
%trials using TSTOOL toolbox
%(http://www.physik3.gwdg.de/tstool/index.html)

%Inputs:
%pS: (point structure) is a structure with fields:
    %X: a pointset of Nt points, D dimensions and Ntr trials in (points x
    %   dimensions x trials) form
    %Nt: number of points (time indexes) of X
    %D: dimensionality of X
    %Ntr: number of trials of X
%Wth: Theiler window in number of points
%treeMetric: distance metric, either 'euclidean' or 'maximum' norm

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
    %treeMetric: distance metric, either 'euclidian' or 'maximum' norm

    
%Permute x so that it is in the (points x trials x dimensions) form
pS.X = permute(pS.X,[1,3,2]);

%Calculate the maximum x point
Xmax = max(pS.X(:));

%Create a Theiler zone of very large numbers to be put between concatenated 
%trials:
theilerZone = 100*Xmax*ones(Wth,pS.Ntr,pS.D);

%Place the points of the Theiler zone inside the pointset
pS.X(end+(1:Wth),:,:) = theilerZone;

%Concatenate trials, i.e. reshape X to take the form 
%((points+theiler window)*Ntr x dimensions):
pS.X = reshape(pS.X,[(pS.Nt+Wth)*pS.Ntr,pS.D]);

%Get rid of the Theiler zone at the end:
pS.X = pS.X(1:end-Wth,:);


%Construct the tree:
pS.kdTree = nn_prepare(pS.X,treeMetric);


%Construct time-to-tree indexes'look-up table:
IND = 1:pS.Nt;
pS.time2treeIND = nan(pS.Nt,pS.Ntr);
for iTr = 1:pS.Ntr;
    pS.time2treeIND(:,iTr) = IND + (iTr-1)*(pS.Nt + Wth) ;
end


%Construct tree-to-time indexes'look-up table:
%Construct a time index vector
tree2timeIND = IND.';
%Add zeros for the Theiler zone
tree2timeIND = [tree2timeIND; zeros(Wth,1)];
%Repeat for Ntr trials
tree2timeIND = repmat(tree2timeIND,[pS.Ntr,1]);
%Get rid of the last Theiler zone at then end:
pS.tree2timeIND = tree2timeIND(1:end-Wth);

pS.Wth=Wth;
pS.treeMetric = treeMetric;