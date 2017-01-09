function MI = DPcalcMI(k,L,Wth,IND,X,Y,XY,trX,trY,trXY)


%This function calculates transinformation or mutual information index, 
%via the estimator in 

%Lindner M, Vicente R, Priesemann V, Wibral M. 2011. 
%TRENTOOL: a Matlab open source toolbox to analyse information flow in 
%time series data with transfer entropy.
%BMC Neurosci. 12:119.

% Gomez-Herrero G, Wu W, Rutanen K, Miguel SC, Pipa G, Vicente R. 2010. 
% Assessing coupling dynamics from an ensemble of time series. 
% arXiv Prepr.

%Inputs:
%k: number of neighbors
%L: number of points in the pointsets (number of rows of X, Y and XY)
%Wth: Theiler window in sample points
%IND: indices of points of interest
%X: pointset of signal x
%Y: pointset of signal y
%XY: pointset of joint space of signals X and Y
%trX: tree of pointset X
%trY: tree of pointset Y
%trXY: tree of pointset XY

%Output:
%MI: estimated mutual information


%Calculate kth nearest neighbor's distance for the joint space xy for each
%point
%TSTOOL
%[index, distance] = nn_search(pointset, atria, query_indices, k,exclude, epsilon)
[~, distXY] = nn_search(XY,trXY,IND,k,Wth,0);


%Calculate number of nearest neighbor's with dist XY for each of the x and
%y signals and for each point

nnX = zeros(L,1);
nnY= zeros(L,1);

%TSTOOL
%[count, neighbors] = range_search(pointset, atria, query_indices, r,
%exclude)
for i=IND;
    [nnTemp, ~] = range_search(X,trX,i,distXY(i,k)-eps,Wth);
    nnX(i) = nnTemp;
    [nnTemp, ~] = range_search(Y,trY,i,distXY(i,k)-eps,Wth);
    nnY(i) = nnTemp;
end


%Estimate TI
MI = psi(k_th)+psi(L)-mean(psi(nnX+1)+psi(nnY+1));
