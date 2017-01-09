function xP = DPpredleastSqOpt(u,Np,nnInd,x,X,k,D,ind_i,ind_i_plus_u)


%This function calculates a prediction based on linear least-squares optimization

% Faes L, Porta A, Nollo G. 2008. 
% Mutual nonlinear prediction as a tool to evaluate coupling strength and 
% directionality in bivariate time series: Comparison among different 
% strategies based on k nearest neighbors. 
% Phys Rev E. 78:1?11.

%http://en.wikipedia.org/wiki/Linear_regression


%Inputs:
%u: prediction horizon in sample points
%Np: number of predicted points
%nnInd: indexes of k nearest neighbors of each point in the X pointset in a
%   matrix of (points,neighbors)
%x: univariate signal
%X: embedded pointset of x
%k: number of nearest neighbors
%D: dimensionality of the embedding
%ind_i: indices of points that predict
%ind_i_plus_u: indices of points to be predicted

%Outputs:
%xP: the predicted points


%Preallocate regressors' matrix
Regressors = zeros(Np,k*D);

%Form the regressors' matrix
%For each predicting point...
for iP=1:Np;
    
    tempN=[];
    %...and for each of its nearest neighbors...
    for iN = 1:k;
        %...get its coordinates to use as regressors..
        tempN = [tempN, X(nnInd(iP,iN),:)];
    end
    %...and put them into the matrix
    Regressors(iP,:) = tempN;
end
clear X;

%Solve the system to find the regression coefficients
RegressionCoefs = lsqlin(Regressors,x(ind_i_plus_u),[],[]);
clear x;

%Finally, calculate the prediction
xP = sum(repmat(RegressionCoefs.',Np,1) .*Regressors,2);

