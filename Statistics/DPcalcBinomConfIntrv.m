function [pl, pu] = DPcalcBinomConfIntrv(p,clim,N,method,n)

%This code calculates confidence intervals of binomial distribution, 
%following methods discussed in:
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
%WRONG:
%    i.e., in the case of equidistant histogram: n = round(p*N);
%    and in any other case (e.g., if kernel smoothing has been applied)
%    n = round(p.*BinSize*N)
%-method: a string out of 
%    'standard': for standard normal approximation (default)
%    'Wilson': for Wilson method
%    'Agresti–Coull': for Agresti–Coull method
%    'ArcSin': for the arcsine method
%    'logit': for the logit method
%    'Clopper–Pearson': for the Clopper–Pearson method

%Outputs:
%-pl: the lower limit of p, a real number in [-oo 1]
%-pu: the upper limit of p, a real number in [0 +oo]

if nargin<4
    method = 'standard';
end

if nargin<5
    n=round(p*N);
end


%The probability not to choose this bin
q=1-p;

%Calculate the (clim/2)th quantile of the standard normal distribution
a = ((100-clim)/2)/100;
k =  norminv(a);

if strcmpi(method,'Wilson')
    k2 = k^2;
    k22 = k2/2;
    Nk2 = N+k2;
    phat = (n+k22)./(Nk2);
    ci = (k./Nk2).*sqrt(N*p.*q+k22/2);
    pl = phat + ci;
    pu = phat - ci;
    
elseif strcmpi(method,'Agresti–Coull')
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
    
elseif strcmpi(method,'ArcSin')
    if isscalar(n)
        n = repmat(n,size(p));
    end
    phat = (n+3/8)/(N+3/4);
    ArcSin = asin(sqrt(phat));
    ci = 1/2*k*sqrt(N);
    pl = sin(ArcSin-ci).^2;
    pu = sin(ArcSin+ci).^2;
    
elseif strcmpi(method,'logit')
    plogit = log(p./q);
    V = N./(n.*(N-n));
    ci = k*sqrt(V);
    plLogit = plogit+ci;
    puLogit = plogit-ci;
    expplLogit = exp(plLogit);
    exppuLogit = exp(puLogit);
    pl = expplLogit./(1+expplLogit);
    pu = exppuLogit./(1+exppuLogit);

elseif strcmpi(method,'Clopper–Pearson')
    if isscalar(n)
        n = repmat(n,size(p));
    end
    pl = betainv(a/2,n,N-n+1);
    pu = betainv(1-a/2,n+1,N-n);
    
else
    ci = k*sqrt(p.*q/N);
    pl = p+ci;
    pu = p-ci;
end


