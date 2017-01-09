function [c, x] = DPcalcCDF1Dens(x,fun,xi,kernel,width)

%Inputs:
%-x: signals in a [points x signals] form
%-fun: a string for one of the functions below
%      -'quantl': quantiles
%      -'interpol': interpolation of empirical cdf
%      -'kernelquantl': quantiles with kernel smoothing
%      -'kernelinterpol': kernel smoothing 
%-xi: bin centers in x or quantiles
%-kernel: a string or a function handle of the kernel to be used
%-width: the bandwidth of the kernel (a positive real number
%               <= max(x)-min(x))

%Outputs
%-c: the cumulative probability distribution (bins x 1)
%-x: bin centers in x space (bins x 1)

x = x(:);
[c, x] = eval([fun,'(x,xi,kernel,width)']);
%c = c.';


function [c, xi] = quantl(x,c,kernel,width)

xi = quantile(x,c);



function [c, xi] = kernelquantl(x,c,kernel,width)

xi = quantile(x,c);

[c, xi] = ksdensity(x,xi,'kernel',kernel,'width',width,'function','cdf');




function [c, xi] = kernelinterpol(x,xi,kernel,width)

[c, xi] = ksdensity(x,xi,'kernel',kernel,'width',width,'function','cdf');



function [c, xi] = interpol(x,xi,kernel,width)

[tempC, tempX] = ecdf(x);
if (tempX(1)==tempX(2))
    tempX(1) = [];
    tempC(1) = [];
end

c = interp1(tempX, tempC,xi);


