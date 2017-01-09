function [pl, pu, p] = DPcalcPvalConfIntrv(p,clim,N,n,correctP,method)

%This code calculates confidence intervals of p-values that result from a
%permutation test with Monte-Carlo randomized simulation, WITHOUT
%replacement, following methods discussed in:
%Brown, L. D., Cai, T. T., & DasGupta, A. (2001).
%Interval Estimation for a Binomial Proportion. 
%Statistical Science, 16(2), 101–133. doi:10.1214/ss/1009213286

%Inputs:
%-p: the observed p-value, a real number in the interval [0 1]
%-clim: the desired confidence level, a real number in the interval [0 100]
%-N: the number of permutation samples, a positive integer
%-n: the number of parmutian samples that led to values more extreme than
%    the observed value of the statistic, a positive integer <=N.

%Optionally:
%-correctP: a flag that signals whether to correct p or not according to the
%           biased estimation of :
%           Smyth, G. K., & Phipson, B. (2010). 
%           Permutation P-values should never be zero: calculating exact  
%           P-values when permutations are randomly drawn. 
%           Statistical Applications in Genetics and Molecular Biology, 9(1),1–12.)
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
%-p: the p-value possibly corrected


%If the number of permutation samples that led to values more extreme than 
%the observed value is given...
if nargin>4
    if correctP
        %...calculate the biased estimate of the true p-value:
        p = (n+1)/(N+1);
    end
end
q=1-p;

if nargin<6 
    method = 'standard';
end

%Calculate the (1 ? ?/2)th quantile of the standard normal distribution
a = ((100-clim)/2)/100;
k =  norminv(a);

if strcmpi(method,'Wilson')
    k2 = k^2;
    k22 = k2/2;
    Nk2 = N+k2;
    phat = (n+k22)/(Nk2);
    ci = (k/Nk2)*sqrt(N*p*q+k22/2);
    pl = phat + ci;
    pu = phat - ci;
    
elseif strcmpi(method,'Agresti–Coull')
    k2 = k^2;
    k22 = k2/2;
    Nk2 = N+k2;
    nk22 = n+k22;
    phat = nk22/Nk2;
    qhat = 1-phat;
    ci = k*sqrt(phat*qhat/N);
    pl = phat+ci;
    pu = phat-ci;
    
elseif strcmpi(method,'ArcSin')
    phat = (n+3/8)/(N+3/4);
    ArcSin = asin(sqrt(phat));
    ci = 1/2*k*sqrt(N);
    pl = sin(ArcSin-ci).^2;
    pu = sin(ArcSin+ci).^2;
    
elseif strcmpi(method,'logit')
    plogit = log(p/q);
    V = N/(n*(N-n));
    ci = k*sqrt(V);
    plLogit = plogit+ci;
    puLogit = plogit-ci;
    expplLogit = exp(plLogit);
    exppuLogit = exp(puLogit);
    pl = expplLogit/(1+expplLogit);
    pu = exppuLogit/(1+exppuLogit);

elseif strcmpi(method,'Clopper–Pearson')
    pl = betainv(a/2,n,N-n+1);
    pu = betainv(1-a/2,n+1,N-n);
    
else
    ci = k*sqrt(p*q/N);
    pl = p+ci;
    pu = p-ci;
end
