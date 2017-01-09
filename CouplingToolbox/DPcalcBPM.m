function [xP, yP] = DPcalcBPM(u,k,D,predFun,prederrFun,L,Wth,IND,x,y,XY,trXY)

%This function calculates a specific point of a bivariate predictability
%map

%Faes L, Porta A, Nollo G. 2008. 
%Mutual nonlinear prediction as a tool to evaluate coupling strength and 
%directionality in bivariate time series: Comparison among different 
%strategies based on k nearest neighbors. 
%Phys Rev E. 78:1?11.


%Inputs:
%u: prediction horizon in sample points
%k: number of neighbors
%D: dimensionality of the embedding
%predFun: a handle to a prediction function
%prederrFun: a handle to a prediction performance function
%L: number of points in the pointsets (number of rows of XY)
%Wth: Theiler window in sample points
%IND: indices of points of interest
%x: x univariate signal
%y: y univariate signal
%XY: pointset of joint space of signals x and y, that can marginally be
%identical to pointsets X and Y for the cases of self- or cross- predictions
%trXY: tree of pointset XY


%Output:
%xP, yP: prediction performances for x and y


IND_i = IND(end-u); %indices of predicting points
IND_i_plus_u =IND(u+1:end); %indices of predicted points

%---------------------Use this if IND~=1:L---------------------------------
% %Find the indices that can be incremented by u, without getting out of
% %the range of the L indices
% INDind = arrayfun( @(IND_i) any(IND==IND_i+u) , IND );
% IND_i = IND(INDind); %indices of predicting points
% IND_i_plus_u = IND(INDind+u); %indices of predicted points
%---------------------Use this if IND~=1:L---------------------------------



%Find their indexes of the k nearest neighbors only of the predicting
%points
%TSTOOL
%[index, distance] = nn_search(pointset, atria, query_indices, k,exclude, epsilon)
[indexXY, ~] = nn_search(XY,trXY,IND_i,k,Wth,0); %distXY


%Perform predictions:
Np = length(IND_i); %number of predicting/predicted points
xP = predFun(u,Np,indexXY,x,XY,k,D,IND_i,IND_i_plus_u); 
yP = predFun(u,Np,indexXY,y,XY,k,D,IND_i,IND_i_plus_u);


%Calculate prediction's performance
xP = prederrFun(x(IND_i_plus_u),xP);
yP = prederrFun(y(IND_i_plus_u),yP);





