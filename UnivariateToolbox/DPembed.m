function [X L tauTotal] = DPembed(x,D,tau,N)

%This function embeds a univariate signal x via time delay embedding

%Inputs:
%x: univariate signal x
%D: embedding dimension of x
%tau: vector of D elements of relative (successive) embedding time lags 
%N: number of points for x and y 

%Outputs:
%X: the embedded pointset
%L: the number of points in the pointset


%If the embedding includes x...
if D>0
    
    tauTotal = cumsum(tau); %the absolute time lags per dimension
    
    L = N-tauTotal(D); 
    
    X=nan(L,D);
    
    %...embed x per dimension
    for iX = 1:D;
        X(:,iX) = x( tauTotal(D-iX+1) + 1 : end-tauTotal(iX));
    end
    
    tauTotal=tauTotal(D);
    
else
    X=[];
    L=0;
    tauTotal=0;
end

