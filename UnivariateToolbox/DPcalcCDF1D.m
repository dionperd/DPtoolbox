function [c, x] = DPcalcCDF1D(x,fun,xi,M,kernel,width)

%Inputs:
%-x: signals in a [points x signals] form
%-fun: a string for one of the functions below
%      -'quantl': quantiles
%      -'interpol': interpolation of empirical cdf
%      -'kernelquantl': quantiles with kernel smoothing
%      -'kernelinterpol': kernel smoothing 
%-xi: bin centers in x or quantiles
%-M: number of signals (columns of x)
%-kernel: a string or a function handle of the kernel to be used
%-width: the bandwidth of the kernel (a positive real number
%               <= max(x)-min(x))

%Outputs
%-c: the cumulative probability distribution (bins x signals)
%-x: bin centers in x space (bins x signals)

[c, x] = eval([fun,'(x,xi,M,kernel,width)']);
%c = c.';        


function [c, xi] = quantl(x,c,M,kernel,width)

xi = quantile(x,c);

c=repmat(c,1,M);


function [c, xi] = kernelquantl(x,c,M,kernel,width)

xi = quantile(x,c);

for iM = 1:M;
    [c(:,iM), xi(:,iM)] = ksdensity(x(:,iM),xi(:,iM),'kernel',kernel,'width',width,'function','cdf');
end



function [c, xi] = kernelinterpol(x,xi,M,kernel,width)

tempXI = xi;
for iM = 1:M;
    [c(:,iM), xi(:,iM)] = ksdensity(x(:,iM),tempXI,'kernel',kernel,'width',width,'function','cdf');
end


function [c, xi] = interpol(x,xi,M,kernel,width)

for iM = 1:M;
    
    [tempC, tempX] = ecdf(x(:,iM));
    if (tempX(1)==tempX(2))
        tempX(1) = [];
        tempC(1) = [];
    end
    
    c(:,iM) = interp1(tempX, tempC,xi);
end

xi =repmat(xi,1,M);
