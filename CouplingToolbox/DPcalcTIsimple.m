function TI = DPcalcTIsimple(k,L,Wth,IND,X,Y,XY,trX,trY,trXY,r)


%This function calculates transinformation or mutual information index, 
%via a simple estimator based on its definition 

% Lambertz M, Vandenhouten R, Grebe R, Langhorst P. 2000. 
% Phase transitions in the common brainstem and related systems investigated 
% by nonstationary time series analysis. 
% J Auton Nerv Syst. 78:141?157.


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


%Calculate the local probability densites as the 
%number of nearest neighbor's within distance r 
%divided by the total number of points
%for each of the x, y and (x,y) signals and for each point

%TSTOOL
%[count, neighbors] = range_search(pointset, atria, query_indices, r,
%exclude)
[pX, ~]  = range_search(X,trX,IND,r,Wth);
[pY, ~]  = range_search(Y,trY,IND,r,Wth);
[pXY, ~] = range_search(XY,trXY,IND,r,Wth);

pX=pX/L;
pY=pY/L;
pXY=pXY/L;

%Estimate TI
TI = mean(log(pXY./(pX.*pY)));
