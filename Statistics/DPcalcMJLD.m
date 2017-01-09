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
    [pHl(:,iSc), pHu(:,iSc)] = DPcalcBinomConfIntrvAgrestiCoull(pH(:,iSc),cfgJLD.clim,thisN);
    
    %alpha stable pdf equidistant
    pAH(:,iSc) = stblpdf(xH(:,iSc),paramsA(1,iSc),paramsA(2,iSc),paramsA(3,iSc),paramsA(4,iSc));
    [pAHl(:,iSc), pAHu(:,iSc)] = DPcalcBinomConfIntrvAgrestiCoull(pAH(:,iSc),cfgJLD.clim,thisN);
    
    %isohistogram pdf
    [pQ(:,iSc), xQ(:,iSc)] = DPcalcStablePDF_1D(thisJL,'isohisto',cfgJLD.dq,Nbins);
    [pQl(:,iSc), pQu(:,iSc)] = DPcalcBinomConfIntrvAgrestiCoull(pQ(:,iSc),cfgJLD.clim,thisN);
    
    %alpha stable pdf isohistogram
    [pAQ(:,iSc)] = stblpdf(xQ(:,iSc),paramsA(1,iSc),paramsA(2,iSc),paramsA(3,iSc),paramsA(4,iSc));
    [pAQl(:,iSc), pAQu(:,iSc)] = DPcalcBinomConfIntrvAgrestiCoull(pAQ(:,iSc),cfgJLD.clim,thisN);
    
    
    %empirical cdf
    [cE(1:thisN,iSc), xi(1:thisN,iSc)] = DPcalcStableCDF_1D(thisJL,'quantiles',Nbins);
    
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
    RMSEcdf(iSc,1)= sqrt( nanmean( ( cE(:,iSc) - cA(:,iSc) ).^2 ) );
end


function [pl, pu] = DPcalcBinomConfIntrvAgrestiCoull(p,clim,N)

%This code calculates confidence intervals of binomial distribution, 
%following the method 'Agresti–Coull' discussed in:
%Brown, L. D., Cai, T. T., & DasGupta, A. (2001).
%Interval Estimation for a Binomial Proportion. 
%Statistical Science, 16(2), 101–133. doi:10.1214/ss/1009213286

%Inputs:
%-p: the probability distribution, a vector of real numbers in the interval [0 1]
%    it is assumed that p has been normalized such that
%    integral_x(p(x))dx=1 (trapz(x,p)=1)
%-clim: the desired confidence level, a real number in the interval [0 100]
%-N: the number of samples, a positive integer


%NOTE!: in all cases N is the number of samples contained in the support of
%the distribution! So, if a specific portion of the quantiles has been used
%in the context of an isohistogram method, the corresponding N has to be
%used be excluding those samples that fall outside the specific support
%space of the distribution.

%Optionally:
%-n: the (hypothetical) number of samples within each bin assuming a
%    binomial distribution
%    default n = round(p*N) in all cases
%   
%Outputs:
%-pl: the lower limit of p, a real number in [-oo 1]
%-pu: the upper limit of p, a real number in [0 +oo]



if nargin<5
    n=round(p*N);
end


%Calculate the (clim/2)th quantile of the standard normal distribution
a = ((100-clim)/2)/100;
k =  norminv(a);

if isscalar(n)
    n = repmat(n,size(p));
end

k2 = k^2;
k22 = k2/2;
Nk2 = N+k2;
nk22 = n+k22;
phat = nk22./Nk2;
qhat = 1-phat;
ci = k*sqrt(phat.*qhat/N);
pl = phat+ci;
pu = phat-ci;




