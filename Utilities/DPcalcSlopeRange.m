function [S, R2] = DPcalcSlopeRange(x,y,N,Np)

%This function finds the slope at every point i of the x-y curve, by fitting
%a line for the [i-Np i+Np] points

%Inputs:
%-x: x coordinate, a vector of real numbers
%-y: y coordinate, a vector of real numbers
%-N: the length of x and y, positive integer
%-Np: a vector of lengths for the segments, vector of positive odd integers
%-th: a vector of different thresholds for slope variability, a vector of
%      real numbers in the interval (0,1)

%Output:
%-S: the slope for each threshold of err

Nps = length(Np);
%Nth = length(th);

%Initialize S and H
S = nan(N,Nps);
R2  = nan(N,Nps);

for iP = 1:Nps;
    
    thisNp = Np(iP);
    
    %For each point of x,y...
    for iP1 = (thisNp-1)/2+1:1:N-(thisNp-1)/2;
        
        %...calculate the index of the included points:
        ind = iP1-(thisNp-1)/2:iP1+(thisNp-1)/2;
        
        thisX = x(ind);
        thisY  = y(ind);
        %...fit a line for those points...
        p=polyfit(thisX, thisY, 1);
        
        %...get the slope of the line
        S(iP1,iP)=p(1);
        
        %Compute the fitted result
        yfit = polyval(p,thisX);
        
        %Compute the residuals
        yresid = thisY - yfit;
        
        %Square the residuals and total them obtain the residual sum of squares:
        SSresid = sum(yresid.^2);
        
        %Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
        SStotal = (thisNp-1) * var(thisY);
        
        %Compute R2:
        R2(iP1,iP) = 1 - SSresid/SStotal;
    end   
    
end

% leg={};
% for iP=Np;
%     leg=[leg; num2str(iP)];
% end
% 
% inds = (Np(1)-1)/2+1:1:N-(Np(1)-1)/2;
% figure;
% subplot(2,1,1)
% title('Slopes')
% plot(exp(x(inds))*0.004,S(inds,:),'o-');
% legend(leg);
%         
% subplot(2,1,2)
% title('R^2')
% plot(exp(x(inds))*0.004,R2(inds,:),'o-');

