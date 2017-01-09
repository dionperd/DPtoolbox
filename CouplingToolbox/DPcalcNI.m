function [Sxy Hxy Nxy Syx Hyx Nyx] = DPcalcNI(k,Wth,IND,X,Y,trX,trY)


%This function calculates nonlinear interdependance indices

%Quiroga RQ, Kraskov A, Kreuz T, Grassberger P. 2002. 
%Performance of different synchronization measures in real data:
%A case study on electroencephalographic signals. 
%Phys Rev E. 65:1?14.

%Inputs:
%k: number of neighbors
%Wth: Theiler window in sample points
%IND: indices of points of interest
%X: pointset of signal x
%Y: pointset of signal y
%trX: tree of pointset X
%trY: tree of pointset Y


%Output:
%H, S, N: indexes of nonlinear interdependance



%TSTOOL
%[index, distance] = nn_search(pointset, atria, query_indices, k,exclude,
%epsilon)


%Calculate euclidean distances of k nearest neighbors for each point of each signal 
%and find their indexes
[indX, distX] = nn_search(X,trX,IND,k,Wth,0);
[indY, distY] = nn_search(Y,trY,IND,k,Wth,0);

%Calculate euclidean distances of k nearest neighbors of the corresponding
%points in time of the other signal
[~, distXY] = nn_search(X,trX,indY,k,Wth,0);
[~, distYX] = nn_search(Y,trY,indX,k,Wth,0);


%Calculate mean euclidean distances of k nearest neighbors for each point 
Rx = mean(distX,2);
Rxy = mean(distXY,2);
Ry = mean(distY,2);
Ryx = mean(distYX,2);


%Estimate NIs

Sxy = mean(Rx./Rxy);
Syx = mean(Ry./Ryx);

Hxy = mean(log(Rx./Rxy));
Hyx = mean(log(Ry./Ryx));

Nxy = mean((Rx-Rxy)./Rx);
Nyx = mean((Ry-Ryx)./Ry);
