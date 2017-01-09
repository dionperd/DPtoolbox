function [pointStat,  multiStat] = DPcalcStat(mVal,surrVal,pointStatMethod,multiStatfun)

%This function calculates a point statistic according to pointStatMethod
%and a multivariate one, according to the function handle multiVarStat


SIZEm = size(mVal); %the size of the measurement values
N = numel(mVal); %total number of the measurement values
SIZEsurr = size(surrVal); %the size of the surrogate values
Nsurr = SIZEsurr(end); %the number of surrogate values per data point

% %Initialize
% pointStat.mval = zeros(N,1);
% pointStat.surrVal = zeros(N,Nsurr);

%Vectorize the statistic
mVal = mVal(:);

%Collapse extra dimensions of surrVal, if any
if length(SIZEsurr)>2
    surrVal  = reshape(surrVal,[N, Nsurr]);
end

%Calculate the mean...
pointStat.mean = reshape(mean(surrVal,2), SIZEm );
%...the standard deviation...
pointStat.std = reshape(std(surrVal,0,2), SIZEm );
%...the most common value...
pointStat.mode = reshape(mode(surrVal,2), SIZEm );
%...the median...
pointStat.median = reshape(median(surrVal,2), SIZEm );
%...the min...
pointStat.min = reshape(min(surrVal,[],2), SIZEm );
%...and the max...
pointStat.max = reshape(max(surrVal,[],2), SIZEm );
%...of the surrogate values

if strcmpi(pointStatMethod,'z')
    
    pointStat.mVal = reshape( ( mVal -pointStat.mean ) ./pointStat.std, SIZEm );
    
    pointStat.surrVal = reshape( ( surrVal - repmat(pointStat.mean,[1 Nsurr]) ) ./repmat(pointStat.std,[1,Nsurr]), SIZEsurr );
    
elseif strcmpi(pointStatMethod,'t')
    
    pointStat.mVal = reshape( sqrt(Nsurr) * ( mVal -pointStat.mean ) ./pointStat.std, SIZEm );
    
    pointStat.surrVal = reshape( sqrt(Nsurr) * ( surrVal - repmat(pointStat.mean,[1 Nsurr]) ) ./repmat(pointStat.std,[1,Nsurr]), SIZEsurr );
    
else
    pointStat.mVal = reshape( mVal, SIZEm );
    
    pointStat.surrVal = reshape( surrVal, SIZEsurr );
end

mVal = reshape( mVal, SIZEm );
surrVal = reshape( surrVal, SIZEsurr );

if ~isempty(multiStatfun) && (SIZEm(1)>1)   
    multiStat = multiStatfun.fun(mVal,surrVal,pointStat,multiStatfun.params);  
else
    multiStat=[];
end