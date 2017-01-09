function rho2 = DPsqCorr(x,y)

%This function calculates the squared correlation between x and y 

% Faes L, Porta A, Nollo G. 2008. Mutual nonlinear prediction as a tool to 
% evaluate coupling strength and directionality in bivariate time series: 
% Comparison among different strategies based on k nearest neighbors. 
% Phys Rev E. 78:1?11.


%rho2 = ( sum(x.Y) )^2 / ( sum(x^2) + sum(y^2) )
%x,y are vectors or matrices in the form (points, signals)

rho2 = sum(x.*y,1).^2./(sum(x.^2,1).*sum(y.*2,1));