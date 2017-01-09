function CCD = DPcalcCCD(r,u,Wth,IND,XY,trXY)


%This function calculates conditional coupling divergence index, 

% Lambertz M, Vandenhouten R, Grebe R, Langhorst P. 2000. 
% Phase transitions in the common brainstem and related systems investigated 
% by nonstationary time series analysis. 
% J Auton Nerv Syst. 78:141?157.


%Inputs:
%r: radius for range search of nearest neighbors
%u: number of time steps ahead for the calculation of the coupling
%   divergence
%l: length of check for conditional coupling divergence in points
%L: number of points in the pointsets (number of rows of X, Y and XY)
%Wth: Theiler window in sample points
%IND: indices of points of interest (increasing monotonously)
%XY: pointset of joint space of signals X and Y
%trXY: tree of pointset XY


%Output
%CCD: estimated coupling divergence index, 


%---------------------Use this if IND~=1:L---------------------------------
% %Find the indices that can be incremented by u, without getting out of
% %the range of the IND indices
% 
% INDind = arrayfun( @(IND_i) any(IND==IND_i+u) , IND );
% IND_i = IND(INDind);
% IND_i_plus_u = IND(INDind+u);
% 
% %Keep only those indices and their increments
% IND = union(IND_i,IND_i_plus_u);
%---------------------Use this if IND~=1:L---------------------------------



%Find the nearest neighbors within distance r 

%TSTOOL
%[count, neighbors] = range_search(pointset, atria, query_indices, r,
%exclude)
%neighbors is a (number of query points x 2)
%the first cell column is a vector of the indices of the nearest neighbors
%the second cell column, gives their distances
[nncountXY, nnXY] = range_search(XY,trXY,IND,r,Wth);


%Calculate the common neighbors between points i and i+u
countXYcommonNN = zeros(length(IND_i),1);

%For each query point i...
for i=1:L-u; %use 'IND_i' instead of '1:L-u' if IND~=1:L;
    %...find the number of common neighbors with the i+u point...
    countXYcommonNN(i) = numel( intersect( nnXY{i,1}+u , nnXY{i+u,1} ) );
end


%Calculate CCD
CCD = mean(countXYcommonNN./nncountXY);


