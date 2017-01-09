function [TExy TEyx] = DPcalcTE(k,IND,L,Wth,X,Y,XY,Xx,Yy,XYx,YXy,trX,trY,trXY,trXx,trYy,trXYx,trYXy)


%This function calculates transfer entropy indices via the estimator in 

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
%X: pointset of signal x
%Y: pointset of signal y
%XY: pointset of joint space of signals X and Y
%trX: tree of pointset X
%trY: tree of pointset Y
%trXY: tree of pointset XY
%Xx: pointset of signal x plus future point of x
%Yy: pointset of signal y plus future point of y
%XYx: pointset of joint space of signals X and Y plus future point of x
%YXy: pointset of joint space of signals X and Y plus future point of y
%trXx: tree of pointset X plus future point of y
%trYy: tree of pointset Y plus future point of y
%trXYx: tree of pointset XY plus future point of x
%trYXy: tree of pointset XY plus future point of y

%Output:
%TExy: x->y estimated transfer entropy
%TEyx: y->x estimated transfer entropy



%Calculate kth nearest neighbor's distance for the joint space xy plus future points of x or y for each
%point

%TSTOOL
%[index, distance] = nn_search(pointset, atria, query_indices, k,exclude, epsilon)
[~, distXYx] = nn_search(XYx,trXYx,IND,k,Wth,0);
[~, distYXy] = nn_search(YXy,trYXy,IND,k,Wth,0);


%Calculate number of nearest neighbor's with dist XY for each of foloowing signals and for each point
nnX = zeros(L,1); %signal x
nnY = zeros(L,1); %signal y
nnXx = zeros(L,1); %signal x plus future point of x
nnYy = zeros(L,1); %signal y plus future point of y
nnXY = zeros(L,1); %joint space of signals x and y for distXYx
nnYX = zeros(L,1); %joint space of signals x and y for distYXy

%TSTOOL
%[count, neighbors] = range_search(pointset, atria, query_indices, r,
%exclude)
for i=IND;
    [nnTemp, ~] = range_search(X,trX,i,distXYx(i,k)-eps,Wth);
    nnX(i) = nnTemp;
    [nnTemp, ~] = range_search(Xx,trXx,i,distXYx(i,k)-eps,Wth);
    nnXx(i) = nnTemp;
    [nnTemp, ~] = range_search(XY,trXY,i,distXYx(i,k)-eps,Wth);
    nnXY(i) = nnTemp;
    [nnTemp, ~] = range_search(Y,trY,i,distYXy(i,k)-eps,Wth);
    nnY(i) = nnTemp;
    [nnTemp, ~] = range_search(Yy,trYy,i,distYXy(i,k)-eps,Wth);
    nnYy(i) = nnTemp;
    [nnTemp, ~] = range_search(XY,trXY,i,distYXy(i,k)-eps,Wth);
    nnYX(i) = nnTemp;
end


%Estimate TE
TExy = psi(k_th)+mean(psi(nnY+1)-psi(nnYy+1)-psi(nnXY+1));
TEyx = psi(k_th)+mean(psi(nnX+1)-psi(nnXx+1)-psi(nnYX+1));