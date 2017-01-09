function [PTExyz PTEyxz] = DPcalcPTE(k,IND,L,Wth,XZ,YZ,XYZ,XZx,YZy,XYZx,YXZy,trXZ,trYZ,trXYZ,trXZx,trYZy,trXYZx,trYXZy)


%This function calculates partial transfer entropy indices via the estimator in 

%Lindner M, Vicente R, Priesemann V, Wibral M. 2011. 
%TRENTOOL: a Matlab open source toolbox to analyse information flow in 
%time series data with transfer entropy.
%BMC Neurosci. 12:119.

% Gomez-Herrero G, Wu W, Rutanen K, Miguel SC, Pipa G, Vicente R. 2010. 
% Assessing coupling dynamics from an ensemble of time series. 
% arXiv Prepr.


%Inputs:
%k: number of neighbors
%L: number of points in the pointsets (number of rows of all pointsets below)
%Wth: Theiler window in sample points
%IND: indices of points of interest
%XZ: pointset of joint signal xz
%YZ: pointset of joint signal yz
%XYZ: pointset of joint space of signals X, Y and Z
%trXZ: tree of pointset XZ
%trYZ: tree of pointset YZ
%trXYZ: tree of pointset XYZ
%XZx: pointset of joint signal xz plus future point of x
%YZy: pointset of joint signal yz plus future point of y
%XYZx: pointset of joint space of signals X, Y and Z plus future point of x
%YXZy: pointset of joint space of signals X, Y and Z plus future point of y
%trXZx: tree of pointset X plus future point of y
%trYZy: tree of pointset Y plus future point of y
%trXYZx: tree of pointset XY plus future point of x
%trYXZy: tree of pointset XY plus future point of y

%Output:
%PTExyz: x->y/z partial transfer entropy
%PTEyxz: y->x/z partial transfer entropy



%Calculate kth nearest neighbor's distance for the joint space xyz plus future points of x or y for each
%point

%TSTOOL
%[index, distance] = nn_search(pointset, atria, query_indices, k,exclude, epsilon)
[~, distXYZx] = nn_search(XYZx,trXYZx,IND,k,Wth,0);
[~, distYXZy] = nn_search(YXZy,trYXZy,IND,k,Wth,0);


%Calculate number of nearest neighbor's with dist XY for each of foloowing signals and for each point
nnXZ = zeros(L,1); %joint signal xz
nnYZ = zeros(L,1); %joint signal yz
nnXZx = zeros(L,1); %joint signal xz plus future point of x
nnYZy = zeros(L,1); %joint signal yz plus future point of y
nnXYZ = zeros(L,1); %joint space of signals x, y and z for distXYZx
nnYXZ = zeros(L,1); %joint space of signals x, y and z for distYXZy

%TSTOOL
%[count, neighbors] = range_search(pointset, atria, query_indices, r,
%exclude)
for i=IND;
    [nnTemp, ~] = range_search(XZ,trXZ,i,distXYZx(i,k)-eps,Wth);
    nnXZ(i) = nnTemp;
    [nnTemp, ~] = range_search(XZx,trXZx,i,distXYZx(i,k)-eps,Wth);
    nnXZx(i) = nnTemp;
    [nnTemp, ~] = range_search(XYZ,trXYZ,i,distXYZx(i,k)-eps,Wth);
    nnXYZ(i) = nnTemp;
    [nnTemp, ~] = range_search(YZ,trYZ,i,distYXZy(i,k)-eps,Wth);
    nnYZ(i) = nnTemp;
    [nnTemp, ~] = range_search(YZy,trYZy,i,distYXZy(i,k)-eps,Wth);
    nnYZy(i) = nnTemp;
    [nnTemp, ~] = range_search(XYZ,trXYZ,i,distYXZy(i,k)-eps,Wth);
    nnYXZ(i) = nnTemp;
end


%Estimate PTE
PTExyz = psi(k_th)+mean(psi(nnYZ+1)-psi(nnYZy+1)-psi(nnXYZ+1));
PTEyxz = psi(k_th)+mean(psi(nnXZ+1)-psi(nnXZx+1)-psi(nnYXZ+1));

