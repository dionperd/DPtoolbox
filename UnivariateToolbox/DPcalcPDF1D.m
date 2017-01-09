function [p, x] = DPcalcPDF1D(x,fun,xi,M,norm,kernel,width)

%Inputs:
%-x: signals in a [points x signals] form
%-fun: a string for one of the functions below
%      -'isodist': equidistant histogram
%      -'isohist': iso-probability histogram
%      -'kernelisodist': equidistant histogram with kernel smoothing
%      -'isohist': iso-probability histogram with kernel smoothing
%-xi: bin centers in x or quantiles
%-M: number of signals (columns of x)
% -norm: equal to 1 in order to normalize pdf to have a total area of
%        1, otherwise equal to 0, default is 1
%-kernel: a string or a function handle of the kernel to be used
%-width: the bandwidth of the kernel (a positive real number
%               <= max(x)-min(x))

%Outputs
%-p: the probability distribution (bins x signals)
%-xi: bin centers in x space (bins x signals)

[p,x] = eval([fun,'(x,xi,M,kernel,width)']);

if norm
    for iM = 1:M;
        p(:,iM) = p(:,iM)./trapz(x(:,iM),p(:,iM));  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab
    end
end


function [p, xi] = isodist(x,xi,M,kernel,width)

p = hist(x,xi);
if size(p,1)==1
    p=p.';
end
%p=p.';

%p = p./repmat(trapz(xi,p),length(xi),1);  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab

xi = repmat(xi,1,M);



function [p, xi] = isohist(x,p,M,kernel,width)

%Quantiles
xi = quantile(x,p);

diffXI = diff(xi);

%The probability distribution as the derivative of the cumulative one...
p = repmat(diff(p),1,M)./diffXI;
%p=p.';

%Bin centers
xi = xi(1:(end-1),:) + diffXI;




function [p, xi] = kernelisodist(x,xi,M,kernel,width)

p = nan(length(xi),M);

for iM = 1:M;
    [p(:,iM), xi(:,iM)] = ksdensity(x(:,iM),xi,'kernel',kernel,'width',width);
end
%p=p.';




function [p, xi] = kernelisohist(x,p,M,kernel,width)

xi = quantile(x,p);

diffXI = diff(xi);

for iM = 1:M;
    [p(:,iM), xi(:,iM)] = ksdensity(x(:,iM),xi(:,iM),'kernel',kernel,'width',width,'function','cdf');
end

%The probability distribution as the derivative of the cumulative one...
p = diff(p)./diffXI;
%p=p.';

%Bin centers
xi = xi(1:(end-1),:) + diffXI;

