function [p, pUnCorr, sign] = DPcalcSignif(stat,alpha,tail,corrMultComp)

%This function decides about the significance of a specific statistic in a
%surrogate test

%Inputs:
%   -stat: a structure containing the statistic, it contains the fields:
%          -mVal: the measurement value array
%          -surrVal: the surrogate values' array
%   -alpha: the alpha value
%   -tail: 1 or 2, for one or two tailed test
%   -corrMultComp: method of correction for multiple comparisons, '','BONF' or
%                  'FDR', default:'' for no correction

%Outputs: 
%   -p: an array of p-values after correction for multiple comparisons
%   -pUnCorr: an array of p-values before correction for multiple
%             comparisons (p=pUnCorr if there is no correction)
%   -sign: an array of 1s and 0s for the significant and non significant
%          data points, respectively

%p, pUnCorr and sign have the same structure as mVal.



SIZEm = size(stat.mVal); %the size of the measurement values
N = numel(stat.mVal); %total number of the measurement values
SIZEsurr = size(stat.surrVal); %the size of the surrogate values
Nsurr = SIZEsurr(end); %the number of surrogate values per data point
%NsurrS = repmat(Nsurr,[N,1]); %minimum p value

%Initialize outputs:
p=nan(N,1);
pUnCorr=nan(N,1);
sign=nan(N,1);

%Vectorize the statistic
stat.mVal = stat.mVal(:);
mVals = repmat(stat.mVal, [1, Nsurr]);

%Collapse extra dimensions of surrVal, if any
if length(SIZEsurr)>2
    stat.surrVal  = reshape(stat.surrVal,[N, Nsurr]);
end

if (tail==1) %For one tailed test...
    
    %...calculate the probability of each data point to be smaller than the
    %surrogate values, with a minimum of 1/Nsurr...
    pUnCorr = max( 1, sum( mVals<stat.surrVal ,2) )/ Nsurr;
    
    
else %For two tailed test...
    
    %...calculate the probability of each data point to be larger or smaller than the
    %surrogate values, with a minimum of 1/Nsurr...
    pUnCorr = max( 1, min( sum( mVals>stat.surrVal ,2),  sum( mVals<stat.surrVal ,2) ) )/Nsurr;
    
    %...divide alpha value by 2 for two-tailed test...
    alpha = alpha/2;
    
end


if strcmpi(corrMultComp,'BONF') %BONF
    
    %...p equals pUnCorr multiplied by number of comparisons...
    p = min(pUnCorr*N, 1);
    
    
elseif strcmpi(corrMultComp,'FDR') % 'FDR'
    
    %...sort p values...
    [pSort, sortIndx] = sort(pUnCorr);
    
    %...FDR correction...
    NN=1./(1:N).'; %a useful constant
    p= sum(NN)*pSort.*(N.*NN); %correction
    p = min(p, 1);
    clear pSort;
    
    %...unsort back...
    [dummy,unsortIndx]=sort(sortIndx);
    p = p(unsortIndx);
    
    
else %...if there is no correction for multiple comparisons...
    
    %...p equals pUnCorr...
    p=pUnCorr;
    
    
end
    
%Calculate significance...
sign = p<=alpha;

%Reshape results back to mVal size:
p= squeeze(reshape(p,SIZEm));
pUnCorr=squeeze(reshape(pUnCorr,SIZEm));
sign=squeeze(reshape(sign,SIZEm));


