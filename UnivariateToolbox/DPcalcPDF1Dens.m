function [p, x] = DPcalcPDF1Dens(x,fun,xi,norm,kernel,width)

%Inputs:
%-x: signals in a [points x signals] form
%-fun: a string for one of the functions below
%      -'isodist': equidistant histogram
%      -'isohist': iso-probability histogram
%      -'kernelisodist': equidistant histogram with kernel smoothing
%      -'isohist': iso-probability histogram with kernel smoothing
%-xi: bin centers in x or quantiles
%-norm: equal to 1 in order to normalize pdf to have a total area of
%       1, otherwise equal to 0, default is 1
%-kernel: a string or a function handle of the kernel to be used
%-width: the bandwidth of the kernel (a positive real number
%               <= max(x)-min(x))

%Outputs
%-p: the probability distribution (bins x 1)
%-xi: bin centers in x space (bins x 1)

x = x(:);
[p,x] = eval([fun,'(x,xi,kernel,width)']);
%p=p.';
if norm
    p = p./trapz(x,p);  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab
end

function [p, xi] = isodist(x,xi,kernel,width)

p = hist(x,xi);
if size(p,1)==1
    p=p.';
end

%p = p./repmat(trapz(xi,p),length(xi),1);  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab







function [p, xi] = isohist(x,p,kernel,width)

%Quantiles
xi = quantile(x,p);

diffXI = diff(xi);

%The probability distribution as the derivative of the cumulative one...
p = diff(p)./diffXI;

%Bin centers
xi = xi(1:(end-1),:) + diffXI;




function [p, xi] = kernelisodist(x,xi,kernel,width)

p = nan(length(xi),1);

p = ksdensity(x,xi,'kernel',kernel,'width',width);




function [p, xi] = kernelisohist(x,p,kernel,width)

xi = quantile(x,p);

diffXI = diff(xi);

[p, xi] = ksdensity(x,xi,'kernel',kernel,'width',width,'function','cdf');


%The probability distribution as the derivative of the cumulative one...
p = diff(p)./diffXI;
         
%Bin centers
xi = xi(1:(end-1),:) + diffXI;

