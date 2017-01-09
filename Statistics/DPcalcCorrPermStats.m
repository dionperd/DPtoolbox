function [C, p, ul, ll, s] = DPcalcCorrPermStats(x,y,CorrType,Nperm,clim,alpha)

[N D] = size(x);

%Calculate correlation for original data and the desired correlation type
%among Pearson', 'Spearman', and 'Kendall'
C = corr(x,y,'type',CorrType);

%Initialize a matrix of Nperm repetitions of the original data
%and arrange it into Nperm positions of a cell
xP = repmat({x},1,Nperm);

%Define a metrix of permuting each column of the original data
funP = @(x) x(randperm(N));

%Apply this function to each column-cell
xP = cellfun(funP,xP,'UniformOutput',0);

%Define a function of correlating each of the permuted data columns with 
%the original y data
funC = @(x,y) corr(x,y,'type',CorrType);

%Create a cell of original y data repetitions
yP = repmat({y},[1,Nperm]);

%Apply the correlation function to all permutations
Cp = cellfun(funC,xP,yP,'UniformOutput',0)';

%Reconvert into a matrix of correlations
Cp = cell2mat(Cp);

%Calculate the lower confidence limit
ll = quantile(Cp,1-clim);

%Calculate the upper confidence limit
ul = quantile(Cp,clim);
      
%Calculate p-value:
plim = 1/Nperm; %lower possible given the number of permutations
pR = sum(Cp>repmat(C,[Nperm, 1]))*plim; %right tail test
pL = sum(Cp<repmat(C,[Nperm, 1]))*plim; %left tail test
p = min(pR,pL); %choose the best one
p(p==0) = plim; 


%Optionally, determine significance:
if (nargout>4) && (nargin>6)
    
    alpha2 = alpha/2; %due to the two tailed test
    
    %if possible given the number of permutations...     
    if alpha2<plim
        s=nan(1,D);
    else
        s = p<=alpha2;
    end
    
else
    s=[];
end

