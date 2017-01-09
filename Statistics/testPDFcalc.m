%The data series
y=randn(10000,1);

%The probability distribution calculation
Nbins = 100;
method = 'isohisto'; 
[p, x, xi] = DPcalcPDF_1D(y,method,Nbins,0);
trapz(x,p)

%The binomial confidence intervals
clim = 95;
N=size(y,1);
%BinSiz = diff(xi);
%WRONG:
%n = round(p.*BinSiz*N); %for isohisto(kernel)
n = round(p*N);%for all cases, WRONG: for equidist(kernel)
method = 'Clopper–Pearson';
[pl, pu] = DPcalcBinomConfIntrv(p,clim,N,n,method);

%Plotting
Color = 'b';
ColorAlpha = 0.25;
figure
hp=patch([x; x(end:-1:1)].',[pu; pl(end:-1:1)].',ones(1,2*length(x)),'edgecolor',Color,'facecolor',Color,'facealpha',ColorAlpha,'edgealpha',0);
hasbehavior(hp,'legend',false);
hold on;
plot(x,p,'color',Color,'LineStyle','-','linewidth',1);
hold off;