function [paramsA, pAQ, pAQl, pAQu, pAH,  pAHl, pAHu,peakAQ, peakAH,...
          pQ, xQ, pQl, pQu, peakQ,...
          pH, xH, pHl, pHu, peakH,...
          cE,cA,  xi, ...
          RMSElogpdfQ, RMSElogpdfH, RMSEcdf] = DPcalcMJLD(JL,cfgJLD)

%Inputs:
%

%Outputs
%
Nt = size(JL,1);
Nsc = cfgJLD.Nsc;
Nbins = cfgJLD.Nbins;

paramsA = nan(4,Nsc);
pAQ = nan(Nbins,Nsc);
pAQl = nan(Nbins,Nsc);
pAQu = nan(Nbins,Nsc);

pAH = nan(Nbins,Nsc);
pAHl = nan(Nbins,Nsc);
pAHu = nan(Nbins,Nsc);

pH = nan(Nbins,Nsc);
pHl = nan(Nbins,Nsc);
pHu = nan(Nbins,Nsc);
xH = nan(Nbins,Nsc);

pQ = nan(Nbins,Nsc);
pQl = nan(Nbins,Nsc);
pQu = nan(Nbins,Nsc);
xQ = nan(Nbins,Nsc);

cA = nan(Nt,Nsc);
xi = nan(Nt,Nsc);
cE = nan(Nt,Nsc);

RMSElogpdfQ = zeros(Nsc,1);
RMSElogpdfH = zeros(Nsc,1);
RMSEcdf = zeros(Nsc,1);

peakAQ = zeros(Nsc,1);
peakAH = zeros(Nsc,1);
peakH = zeros(Nsc,1);
peakQ = zeros(Nsc,1);


for iSc = 1:Nsc; %...and for each scale...
    
    %Get
    thisJL = sort( JL(~isnan(JL(:,iSc)) ,iSc) );
    thisN = length(thisJL);
    
    %alpha stable fitting
    paramsA(:,iSc) = stblfit(thisJL);   
    
    %equidistant histogram pdf
    [pH(:,iSc), xH(:,iSc)] = DPcalcStablePDF_1D(thisJL,'equidist',cfgJLD.dq,Nbins);
    [pHl(:,iSc), pHu(:,iSc)] = DPcalcBinomConfIntrv(pH(:,iSc),cfgJLD.clim,thisN,'AgrestiCoull');
    
    %alpha stable pdf equidistant
    pAH(:,iSc) = stblpdf(xH(:,iSc),paramsA(1,iSc),paramsA(2,iSc),paramsA(3,iSc),paramsA(4,iSc));
    [pAHl(:,iSc), pAHu(:,iSc)] = DPcalcBinomConfIntrv(pAH(:,iSc),cfgJLD.clim,thisN,'AgrestiCoull');
    
    %isohistogram pdf
    [pQ(:,iSc), xQ(:,iSc)] = DPcalcStablePDF_1D(thisJL,'isohisto',cfgJLD.dq,Nbins);
    [pQl(:,iSc), pQu(:,iSc)] = DPcalcBinomConfIntrv(pQ(:,iSc),cfgJLD.clim,thisN,'AgrestiCoull');
    
    %alpha stable pdf isohistogram
    [pAQ(:,iSc)] = stblpdf(xQ(:,iSc),paramsA(1,iSc),paramsA(2,iSc),paramsA(3,iSc),paramsA(4,iSc));
    [pAQl(:,iSc), pAQu(:,iSc)] = DPcalcBinomConfIntrv(pAQ(:,iSc),cfgJLD.clim,thisN,'AgrestiCoull');
    
    
    %empirical cdf
    [cE(1:thisN,iSc), xi(1:thisN,iSc)] = DPcalcStableCDF_1D(thisJL,'empirical');
    
    %alpha stable cdf
    cA(1:thisN,iSc) = stblcdf(xi(1:thisN,iSc),paramsA(1,iSc),paramsA(2,iSc),paramsA(3,iSc),paramsA(4,iSc));    
    
    
    % and the corresponding pdf peaks
    
    [dummy, ind] = max(pAQ(:,iSc));
    peakAQ(iSc,1) = xQ(ind,iSc);
    
    [dummy, ind] = max(pAH(:,iSc));
    peakAH(iSc,1) = xH(ind,iSc);
    
    [dummy, ind] = max(pH(:,iSc));
    peakH(iSc,1) = xH(ind,iSc);
    
    [dummy, ind] = max(pQ(:,iSc));
    peakQ(iSc,1) = xQ(ind,iSc);
    
    
    %RMSE of non-parametric - parametric isohistogram pdf
    ind = (pQ(:,iSc)>0) & (pQ(:,iSc)<inf) & (pAQ(:,iSc)>0) & (pAQ(:,iSc)<inf);
    RMSElogpdfQ(iSc,1) = sqrt( nanmean( ( log(pQ(ind,iSc)) - log(pAQ(ind,iSc)) ).^2 ) );
    
    %RMSE of non-parametric - parametric equidistant pdf
    ind = (pH(:,iSc)>0) & (pH(:,iSc)<inf) & (pAH(:,iSc)>0) & (pAH(:,iSc)<inf);
    RMSElogpdfH(iSc,1) = sqrt( nanmean( ( log(pH(ind,iSc)) - log(pAH(ind,iSc)) ).^2 ) );
    
    %RMSE of empirical - parametric cdf
    ind = (cE(:,iSc)>=0) & (cE(:,iSc)<=1) & (cA(:,iSc)>=0) & (cA(:,iSc)<=1);
    RMSEcdf(iSc,1)= sqrt( nanmean( ( cE(ind,iSc) - cA(ind,iSc) ).^2 ) );
end




