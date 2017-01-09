
%This function is a naive calculator of 1D Shannon entropy

%-------------------Very simple code made by DP----------------------------
%function NE = DPcalcNE(x,Nbins)
% %Calculate a !D histogram of Nbins bins and devide by length(x) to get the
% %probability distribution:
% x=x(:);
% N=length(x);
% p = hist(x,Nbins)/N;
% 
% %Calculate entropy normalized with the number of bins, in order to give
% %values between 0 and 1:
% NE = - sum(p(p~=0).*log2(p(p~=0)))/log2(Nbins);
%-------------------Very simple code made by DP----------------------------



function [NE, estimation, nbias,sigma]=DPcalcNEmulti(x,ncell,approach,base,normInp,normOut)
%ENTROPY   Estimates the entropy of stationary signals with
%          independent samples using various approaches.
%   [ESTIMATION,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X) or
%   [ESTIMATION,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR) or
%   [ESTIMATION,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR,APPROACH) or
%   [ESTIMATION,NBIAS,SIGMA,DESCRIPTOR] = ENTROPY(X,DESCRIPTOR,APPROACH,BASE)
%
%   ESTIMATION    : The entropy estimate
%   NBIAS       : The N-bias of the estimate
%   SIGMA       : The standard error of the estimate
%   DESCRIPTOR  : The descriptor of the histogram, seel alse ENTROPY
%
%   X           : The time series to be analyzed, a row vector
%   DESCRIPTOR  : Where DESCRIPTOR=[LOWERBOUND,UPPERBOUND,NCELL]
%     LOWERBOUND: Lowerbound of the histogram
%     UPPERBOUND: Upperbound of the histogram
%     NCELL     : The number of cells of the histogram       
%   APPROACH    : The method used, one of the following ones:
%     'unbiased': The unbiased estimate (default)
%     'mmse'    : The minimum mean square error estimate
%     'biased'  : The biased estimate
%   BASE        : The base of the logarithm; default e
%
%   See also: http://www.cs.rug.nl/~rudy/matlab/

%   R. Moddemeijer 
%   Copyright (c) by R. Moddemeijer
%   $Revision: 1.1 $  $Date: 2001/02/05 08:59:36 $

% Moddemeijer, R. 
% On Estimation of Entropy and Mutual Information of Continuous Distributions, 
% Signal Processing, 1989, vol. 16, nr. 3, pp. 233-246, abstract , BibTeX ,
% 
% For the principle of Minimum Mean Square Error estimation see:
% Moddemeijer, R. 
% An efficient algorithm for selecting optimal configurations of AR-coefficients, 
% Twentieth Symp. on Information Theory in the Benelux, May 27-28, 1999, Haasrode (B), pp 189-196, 
% eds. A. Barbé et. al., Werkgemeenschap Informatie- en Communicatietheorie, Enschede (NL), 
% and IEEE Benelux Chapter on Information Theory, ISBN: 90-71048-14-4, abstract , BibTeX ,


%-----Modified by DP, look in external folder for the original code--------


%Normalize:
if normInp
   x=zscore(x);
end

[N, D] = size(x);

%%Vectorize
%x=x(:);

%histogram parameters
%set the limits of the histogram as [-3*std 3*std]
lowerbound=-3;
upperbound=3;
range = 6;
logRange = log(range);

%Calculate histogram
h=histogram(x,lowerbound, upperbound, ncell,N,D); 




% %Initialize
% estimation=0; %entropy estimate
% sigma=0;
% count=0; %number of points

%Calculate probability distribution
% for n=1:ncell
%   if h(n)~=0 
%     logf=log(h(n));
%   else
%     logf=0;
%   end;
%   count=count+h(n);
%   estimation=estimation-h(n)*logf;
%   sigma=sigma+h(n)*logf^2;
% end;

%DP faster code:
h0 = h(h~=0);
count = sum(h0); %the total number of points in the calculation
estimation = -sum( h0.*log(h0) ); 
sigma = sum( h0.*log(h0).^2 );


% biased estimate
NbinsBias = ncell*D;
estimation=estimation/count;
sigma   =sqrt( (sigma/count-estimation^2)/(count-1) );
estimation=estimation+log(count)+logRange- log(NbinsBias);
nbias   =-(NbinsBias-1)/(2*count);

% conversion to unbiased estimate
if approach(1)=='u'
  estimation=estimation-nbias;
  nbias=0;
end;

% conversion to minimum mse estimate
if approach(1)=='m'
  estimation=estimation-nbias;
  %nbias=0; DP: unnecessary
  lambda=estimation^2/(estimation^2+sigma^2);
  nbias   =(1-lambda)*estimation;
  estimation=lambda*estimation;
  sigma   =lambda*sigma;
end;

% base transformation
logBase = log(base);
estimation=estimation/logBase;
nbias   =nbias   /logBase;
sigma   =sigma   /logBase;

if normOut
    %DP:transform estimation to the range [0 1]
    NE = estimation/(logRange/log(base));
end

function result=histogram(x,lower, upper, ncell,N,D)
%HISTOGRAM Computes the frequency histogram of the row vector x.
%   [RESULT,DESCRIPTOR] = HISTOGRAM(X) or
%   [RESULT,DESCRIPTOR] = HISTOGRAM(X,DESCRIPTOR) or
%where
%   DESCRIPTOR = [LOWER,UPPER,NCELL]
%
%   RESULT    : A row vector containing the histogram
%   DESCRIPTOR: The used descriptor
%
%   X         : The row vector be analyzed
%   DESCRIPTOR: The descriptor of the histogram
%     LOWER   : The lowerbound of the histogram
%     UPPER   : The upperbound of the histogram
%     NCELL   : The number of cells of the histogram
%
%   See also: http://www.cs.rug.nl/~rudy/matlab/

%   R. Moddemeijer 
%   Copyright (c) by R. Moddemeijer
%   $Revision: 1.2 $  $Date: 2001/02/05 09:54:29 $
%-----Modified by DP, look in external folder for the original code--------


%DP: this part should be transferred in the parameter check of the cfg
%preparatory function:
% if nargin==1
%   minx=min(x);
%   maxx=max(x);
%   delta=(maxx-minx)/(length(x)-1);
%   ncell=ceil(sqrt(length(x)));
%   descriptor=[minx-delta/2,maxx+delta/2,ncell];
% end;

% lower=descriptor(1);
% upper=descriptor(2);
% ncell=descriptor(3);

% if ncell<1 
%   error('Invalid number of cells')
% end;
% 
% if upper<=lower
%   error('Invalid bounds')
% end;

result=zeros(ncell,D);

y=repmat(round( (x-lower)/(upper-lower)*ncell + 1/2 ),1,D);
for iD=1:D;
    for iP=1:N;
        index=y(iP,iD);
        if index >= 1 && index<=ncell
            result(index,iD)=result(index,iD)+1;
        end;
    end;
end
result=result(:);

