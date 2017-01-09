function multiStat = DPcrossCorrStim(mVal,surrVal,pointStat,params)

%This function calculates a multivariate surrogate statistic as the maximum
%cross correlation between either the original time series of the measure 
%and of its surrogates (mVal, surrVal) or of some function of their point 
%statistics (pointStat.mVal, pointStat.surrVal)


SIZEsurr=size(surrVal); %size of surrogate data
D=length(SIZEsurr)-1; %dimension of measurement
Nsurr = SIZEsurr(D+1); %he number of surrogates
if (D==1) %if the statistic is univariate

    %Select the exact point statistic to be used for the multivariate one
    if strcmpi(params.mode, 'max')
        
        stat = mVal - pointStat.max;
        
        statSurr = surrVal - repmat(pointStat.max,[1 Nsurr]);

        
    elseif strcmpi(params.mode, 'mean')
        
        stat = mVal - pointStat.mean;
        
        statSurr = surrVal - repmat(pointStat.mean,[1 Nsurr]);
        
        
    elseif strcmpi(params.mode, 'median')
        
        stat = mVal - pointStat.median;
        
        statSurr = surrVal - repmat(pointStat.median,[1 Nsurr]);
        
        
    elseif strcmpi(params.mode, 'mode')
        
        stat = mVal - pointStat.mode;
        
        statSurr = surrVal - repmat(pointStat.mode,[1 Nsurr]);
        
        
    elseif strcmpi(params.mode, 'pointStat')
        
        stat = pointStat.mVal;
        
        statSurr = pointStat.surrVal;
        
        
    else
        
        stat = mVal;
        
        statSurr = surrVal;
        
        
    end
    
    clear mVal surrVal poinStat; %we don't need these any more
    
    
    %...calculate multivariate statistic for the measure: cross correlation...
    r=xcorr(stat,params.Stim,params.lag,params.bias);
    
    %...select only positive or negative lags if required...
    if strcmpi(params.lagDir,'p')
        %...get rid of negative lags...
        r = r(params.lag+1:end);
    elseif strcmpi(params.lagDir,'n')
        %...get rid of positive lags...
        r = r(params.lag+1:end);
    end
    
    %...find the maximum cross-correlation...
    multiStat.mVal=max(r);
    
    
    %...for each surrogate...
    for iST=1:Nsurr;
        
        %...calculate multivariate statistic: cross correlation...
        r=xcorr(statSurr(:,iST),params.Stim,params.lag,params.bias);
        
        %...select only positive or negative lags if required...
        if strcmpi(params.lagDir,'p')
            %...get rid of negative lags...
            r = r(params.lag+1:end);
        elseif strcmpi(params.lagDir,'n')
            %...get rid of positive lags...
            r = r(params.lag+1:end);
        end
        
        %...find the maximum cross-correlation...
        multiStat.surrVal(iST) = max(r);
        
    end
    
    %Calculate the mean...
    multiStat.mean = mean(multiStat.surrVal);
    %...the standard deviation...
    multiStat.std = std(multiStat.surrVal);
    %...the most common value...
    multiStat.mode = mode(multiStat.surrVal);
    %...the median...
    multiStat.median = median(multiStat.surrVal);
    %...the min...
    multiStat.min = min(multiStat.surrVal);
    %...and the max...
    multiStat.max = max(multiStat.surrVal);
    %...of the multivariate statistic of the surrogate values


    %Calculate the distance of the measures's value to the mean...
    multiStat.meanDist = multiStat.mVal-multiStat.mean;
    %...the most common value...
    multiStat.modeDist = multiStat.mVal-multiStat.mode;
    %...the median...
    multiStat.medianDist = multiStat.mVal-multiStat.median;
    %...the min...
    multiStat.minDist = multiStat.mVal-multiStat.min;
    %...and the max...
    multiStat.maxDist = multiStat.mVal-multiStat.max;
    %...of the multivariate statistic of the surrogate values
    
    %Calculate the normalized distance of the measures's value to the mean...
    multiStat.meanDistNorm = multiStat.meanDist/multiStat.std;
    %...the most common value...
    multiStat.modeDistNorm = multiStat.modeDist/multiStat.std;
    %...the median...
    multiStat.medianDistNorm = multiStat.medianDist/multiStat.std;
    %...the min...
    multiStat.minDistNorm = multiStat.minDist/multiStat.std;
    %...and the max...
    multiStat.maxDistNorm =multiStat.maxDist/multiStat.std;
    %...of the multivariate statistic of the surrogate values
    
    
elseif (D==2) %if the statistic is bivariate
    
        
    %Select the exact point statistic to be used for the multivariate one
    if strcmpi(params.mode, 'max')
        
        stat = mVal - pointStat.max;
        
        statSurr = surrVal - repmat(pointStat.max,[1 1 Nsurr]);
        
        
    elseif strcmpi(params.mode, 'mean')
        
        stat = mVal - pointStat.mean;
        
        statSurr = surrVal - repmat(pointStat.mean,[1 1 Nsurr]);
        
        
    elseif strcmpi(params.mode, 'median')
        
        stat = mVal - pointStat.median;
        
        statSurr = surrVal - repmat(pointStat.median,[1 1 Nsurr]);
        
        
    elseif strcmpi(params.mode, 'mode')
        
        stat = mVal - pointStat.mode;
        
        statSurr = surrVal - repmat(pointStat.mode,[1 1 Nsurr]);
        
        
    elseif strcmpi(params.mode, 'pointStat')
        
        stat = pointStat.mVal;
        
        statSurr = pointStat.surrVal;
        
        
    else
        
        stat = mVal;
        
        statSurr = surrVal;
        
        
    end
    
    clear mVal surrVal poinStat; %we don't need these any more
    
    
    for iT = 1:SIZEsurr(2);
        
        %...calculate multivariate statistic for the measure: cross correlation...
        r=xcorr(stat(:,iT),params.Stim,params.lag,params.bias);
        
        if strcmpi(params.lagDir,'p')
            %...get rid of negative lags...
            r = r(params.lag+1:end);
        elseif strcmpi(params.lagDir,'n')
            %...get rid of negative lags...
            r = r(params.lag+1:end);
        end
        
        %...find the maximum cross-correlation...
        multiStat.mVal(iT)=max(r);
        
        
        %...for each surrogate...
        for iST=1:Nsurr;
            
            %...calculate multivariate statistic: cross correlation...
            r=xcorr(statSurr(:,iT,iST),params.Stim,params.lag,params.bias);
            %...get rid of negative lags...
            r = r(params.lag+1:end);
            
            %...find the maximum cross-correlation...
            multiStat.surrVal(iT,iST) = max(r);
            
        end  %end of iST (surrogate data) loop
        
%         %Calculate the mean...
%         multiStat.meanPtrial(iT) = mean(multiStat.surrVal(iT,:));
%         %...the standard deviation...
%         multiStat.stdPtrial(iT) = std(multiStat.surrVal(iT,:));
%         %...the most common value...
%         multiStat.modePtrial(iT) = mode(multiStat.surrVal(iT,:));
%         %...the median...
%         multiStat.medianPtrial(iT) = median(multiStat.surrVal(iT,:));
%         %...the min...
%         multiStat.minPtrial(iT) = min(multiStat.surrVal(iT,:));
%         %...and the max...
%         multiStat.maxPtrial(iT) = max(multiStat.surrVal(iT,:));
%         %...of the multivariate statistic of the surrogate values per Trial
        
        
    end %end of iT (trial) loop
    
    
    %Calculate the mean...
    multiStat.mean = mean(multiStat.surrVal(:));
    %...the standard deviation...
    multiStat.std = std(multiStat.surrVal(:));
    %...the most common value...
    multiStat.mode = mode(multiStat.surrVal(:));
    %...the median...
    multiStat.median = median(multiStat.surrVal(:));
    %...the min...
    multiStat.min = min(multiStat.surrVal(:));
    %...and the max...
    multiStat.max = max(multiStat.surrVal(:));
    %...of the multivariate statistic of the surrogate values across trials
    
    %Calculate the distance of the measures's value to the mean...
    multiStat.meanDist = multiStat.mVal-multiStat.mean;
    %...the most common value...
    multiStat.modeDist = multiStat.mVal-multiStat.mode;
    %...the median...
    multiStat.medianDist = multiStat.mVal-multiStat.median;
    %...the min...
    multiStat.minDist = multiStat.mVal-multiStat.min;
    %...and the max...
    multiStat.maxDist = multiStat.mVal-multiStat.max;
    %...of the multivariate statistic of the surrogate values
    
    %Calculate the normalized distance of the measures's value to the mean...
    multiStat.meanDistNorm = multiStat.meanDist/multiStat.std;
    %...the most common value...
    multiStat.modeDistNorm = multiStat.modeDist/multiStat.std;
    %...the median...
    multiStat.medianDistNorm = multiStat.medianDist/multiStat.std;
    %...the min...
    multiStat.minDistNorm = multiStat.minDist/multiStat.std;
    %...and the max...
    multiStat.maxDistNorm =multiStat.maxDist/multiStat.std;
    %...of the multivariate statistic of the surrogate values
    
else
    
    multiStat=[];
    
end %if D==1 or D==2



