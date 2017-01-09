function [CPxy CPyx MPxy MPyx PIxy PIyx NPIxy, NPIyx] = DPcalcMP(tau,n,m,treeMetric,u,k,predFun,prederrFun,Wth,IND,x,y)


%This function calculates mutual predictions indices

%Faes L, Porta A, Nollo G. 2008. 
%Mutual nonlinear prediction as a tool to evaluate coupling strength and 
%directionality in bivariate time series: Comparison among different 
%strategies based on k nearest neighbors. 
%Phys Rev E. 78:1?11.


%Inputs:
%tau: embedding time lag
%n,m: vectors of embedding dimensionalities for signals x and y
%   respectively, to be used for the creation of the Bivariate Predictability
%   Map (starting from dimension 1)
%treeMetric: the metric to use for the construction of the tree pointset, 
%   either 'euclidean', or 'maximum'
%u: prediction horizon in sample points
%k: number of neighbors
%predFun: a handle to a prediction function
%prederrFun: a handle to a prediction performance function
%L: number of points in the pointsets (number of rows of XY)
%Wth: Theiler window in sample points
%IND: indices of points of interest
%x: univariate signal x
%y: pointset of signal y



%Outputs:
%CPxy, CPyx: estimated cross-predictions
%MPxy, MPyx: estimated mixed-predictions
%PIxy, PIyx: estimated predictability-improvements
%NPIxy, NPIyx: estimated normalized predictability-improvements


N = numel(x);

nN = numel(n);
mN = numel(m);

%Allocate memory for the Bivariate Predictability Maps
BPMx = zeros(nN,mN);
BPMy = zeros(mN,nN);


%First, calculate for x signal ONLY...
for iN=n;
    
    %..embed x...
    tauX = [0 repmat(tau,1,n(iN)-1)];
    [X, Lx, tauX] = DPembed(x,n(iN),tauX,N);
    
    %...and generate the corresponding tree pointset...
    %TSTOOL
    %atria = nn_prepare(pointset, metric, clustersize)
    trX = nn_prepare(X, treeMetric);
    
    %...and caculate the x self-predictability...
    [BPMx(iN,1), BPMy(1,iN)] = DPcalcBPM(u,k,n(iN),predFun,prederrFun,Lx,Wth,IND,x(tauX+1:end),y(tauX+1:end),X,trX);
end
clear X trX;


%...second, calculate for y signal ONLY...
for iM=m;
    
    %..embed x...
    tauY = [0 repmat(tau,1,m(iM)-1)];
    [Y, Ly, tauY] = DPembed(y,m(iM),tauY,N);
    
    %...and generate the corresponding tree pointset...
    %TSTOOL
    %atria = nn_prepare(pointset, metric, clustersize)
    trY = nn_prepare(Y, treeMetric);
    
    %...and caculate the y self-predictability...
    [BPMx(1,iM+1), BPMy(iM+1,1)] = DPcalcBPM(u,k,m(iM),predFun,prederrFun,Ly,Wth,IND,x(tauY+1:end),y(tauY+1:end),Y,trY);
end
clear Y trY;


%...and finally, for each combination of embedding dimensionalities...
for iN=n;
    
    for iM=m;
        
        %..embed x...
        tauX = [0 repmat(tau,1,n-1)];
        [X, Lx, tauX] = DPembed(x,n,tauX,N);
        
        %..embed y...
        tauY = [0 repmat(tau,1,m-1)];
        [Y, Ly, tauY] = DPembed(y,m,tauY,N);
        
        %...join the pointsets...
        L=min(Lx,Ly);
        tauXY = max(tauX, tauY);
        XY=[X(1:L) Y(1:L)];
        
        clear X Y;
        
        %...and generate the corresponding tree pointset...
        %TSTOOL
        %atria = nn_prepare(pointset, metric, clustersize)
        trXY = nn_prepare(XY, treeMetric);
        
        %...and caculate the predictability...
        [BPMx(iN+1,iM+1), BPMy(iM+1,iN+1)] = DPcalcBPM(u,k,n(iN)+m(iM),predFun,prederrFun,L,Wth,IND,x(tauXY+1:end),y(tauXY+1:end),XY,trXY);
    end
end
clear x y XY trXY IND;


%Calculate cross-predictions:
CPxy = max(BPMx(1,:));
CPyx = max(BPMy(1,:));

%Calculate self-predictions:
SPx = max(BPMx(:,1));
SPy = max(BPMy(:,1));

%Calculate mixed-predictions:
MPxy = max(BPMx(:));
MPyx = max(BPMy(:));

%Calculate predictability-improvements:
PIxy = MPxy - SPx;
PIyx = MPyx - SPy;

%Calculate normalized predictability-improvements:
NPIxy = PIxy/(1 - SPx);
NPIyx = PIyx/(1 - SPy);


