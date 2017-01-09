function [r, p, pol, R2] = DPlinRegress(x,y)

%Compute the correlation coefficient and the p-value
[r, p]= corrcoef(x,y);
r=r(1,2);
p=p(1,2);

%Compute the regression line
[pol, S]= polyfit(x,y,1);

%Compute the fitted result
yfit = polyval(pol,x);

%Compute the residuals
yresid = y - yfit;

%Square the residuals and sum them to obtain the residual sum of squares:
SSresid = sum(yresid.^2);

%Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(y)-1) * var(y);

%Compute R2:
R2 = 1 - SSresid/SStotal;


