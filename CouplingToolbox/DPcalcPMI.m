function PMI = DPcalcPMI(k,L,Wth,IND,XZ,YZ,Z,XYZ,trXZ,trYZ,trZ,trXYZ)


%This function calculates partial transinformation 
%or mutual information index, via the estimator in 

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
%XZ: pointset of joint signal xz
%YZ: pointset of joint signal yz
%Z: pointset of signal z
%XYZ: pointset of joint space of signals X and Y
%trXZ: tree of pointset XZ
%trYZ: tree of pointset YZ
%trZ: tree of pointset Z
%trXYZ: tree of pointset XY

%Output:
%PMI: estimated partial mutual information


%Calculate kth nearest neighbor's distance for the joint space xyz for each
%point
%TSTOOL
%[index, distance] = nn_search(pointset, atria, query_indices, k,exclude, epsilon)
[~, distXYZ] = nn_search(XYZ,trXYZ,IND,k,Wth,0);


%Calculate number of nearest neighbor's with dist XY for each of the z, xz and
%yz signals and for each point

nnXZ = zeros(L,1);
nnYZ= zeros(L,1);
nnZ= zeros(L,1);

%TSTOOL
%[count, neighbors] = range_search(pointset, atria, query_indices, r,
%exclude)
for i=IND;
    [nnTemp, ~] = range_search(XZ,trXZ,i,distXYZ(i,k)-eps,Wth);
    nnXZ(i) = nnTemp;
    [nnTemp, ~] = range_search(YZ,trYZ,i,distXYZ(i,k)-eps,Wth);
    nnYZ(i) = nnTemp;
    [nnTemp, ~] = range_search(Z,trZ,i,distXYZ(i,k)-eps,Wth);
    nnZ(i) = nnTemp;
end


%Estimate CTI
PMI = psi(k_th)-mean(psi(nnXZ+1)+psi(nnYZ+1)-psi(nnZ+1));
