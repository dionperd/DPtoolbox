%testDPcalcStablePDF_1D

%clear all

%Create distribution:

M=1;

%--------------------------------------------------------------------------
%x=randn(N,M);
%--------------------------------------------------------------------------
% x = y(:,31,31);
% sc = 300;
% x = x( (1+sc):end) - x( 1:(end-sc)) ;
%--------------------------------------------------------------------------
%parameters
alpha = 2;
beta = 1;
gamma = 0.03;
delta = 0;

%support
x = [-10:0.1:10].';
indP = x>0;
xlog = log(x(indP));

%number of samples
N = 2500;

%alpha stable:
tic
x = stblrnd(alpha,beta,gamma,delta,[N,1]);
%--------------------------------------------------------------------------

x = unique(x);

N=length(x);

%%Estimate distributions

Nbins = round(sqrt(N));
%Nbins = 100;

%Quantile:
q0  = 1/Nbins


%parametric estimation

%alpha stable:
tic
[pA, xA, paramsA] = DPcalcStablePDF_1D(x,'ASDq',q0);
cA = stblcdf(x,paramsA(1),paramsA(2),paramsA(3),paramsA(4));
%pA = interp1(xA,pA,x);
%xA = x;
toc
paramsA
temp = xA>0;
xAlog = log(xA(temp));
pAlog = log(pA(temp));
[maxP, indP] = max(pA);
xAp = xA(indP);
xAp
[maxP, indP] = max(pAlog);
pAslopes = diff(pAlog(indP:end))./diff(xAlog(indP:end));
inds  = ~isnan(pAslopes) & ~isinf(pAslopes);
pAslopesM = mean(pAslopes(inds))
pAslopesV = std(pAslopes(inds))
pAslopesVM = pAslopesV/abs(pAslopesM)

% %LSD:
% tic
% [pL, xL, paramsL] = DPcalcStablePDF_1D(x,'LSDh',q0);
% cL = levystblcdf(x,paramsL(1),paramsL(2),paramsL(3),paramsL(4));
% pL = interp1(xL,pL,x);
% xL = x;
% toc
% paramsL
% temp = xL>0;
% xLlog = log(xL(temp));
% pLlog = log(pL(temp));
% [maxP, indP] = max(pL);
% xLp = xL(indP);
% xLp
% [maxP, indP] = max(pLlog);
% pLslopes = diff(pLlog(indP:end))./diff(xLlog(indP:end));
% inds  = ~isnan(pLslopes) & ~isinf(pLslopes);
% pLslopesM = mean(pLslopes(inds))
% pLslopesV = std(pLslopes(inds))
% pLslopesVM = pLslopesV/abs(pLslopesM)


%non-parametric estimation

%empirical cumulative distribution:
[cE xE] = ecdf(x);
pE = diff(cE)./diff(xE);
xE = xE(2:end);
temp = xE>0;
xElog = log(xE(temp));
pElog = log(pE(temp));
[maxP, indP] = max(pE);
xEp = xE(indP);
xEp
[maxP, indP] = max(pElog);
pEslopes = diff(pElog(indP:end))./diff(xElog(indP:end));
inds  = ~isnan(pEslopes) & ~isinf(pEslopes);
pEslopesM = mean(pEslopes(inds))
pEslopesV = std(pEslopes(inds))
pEslopesVM = pEslopesV/abs(pEslopesM)

%equidistang histogram calculation
tic
[pH, xH] = DPcalcStablePDF_1D(x,'equidist',q0);
pH = interp1(xH,pH,x);
xH = x;
toc
temp = xH>0;
xHlog = log(xH(temp));
pHlog = log(pH(temp));
[maxP, indP] = max(pH);
xHp = xH(indP);
xHp
[maxP, indP] = max(pHlog);
pHslopes = diff(pHlog(indP:end))./diff(xHlog(indP:end));
inds  = ~isnan(pHslopes) & ~isinf(pHslopes);
pHslopesM = mean(pHslopes(inds))
pHslopesV = std(pHslopes(inds))
pHslopesVM = pHslopesV/abs(pHslopesM)

%iso-histogram calculation
tic
[pQ, xQ] = DPcalcStablePDF_1D(x,'isohisto',q0);
pQ = interp1(xQ,pQ,x);
xQ = x;
toc
temp = xQ>0;
xQlog = log(xQ(temp));
pQlog = log(pQ(temp));
[maxP, indP] = max(pQ);
xQp = xQ(indP);
xQp
[maxP, indP] = max(pQlog);
pQslopes = diff(pQlog(indP:end))./diff(xQlog(indP:end));
inds  = ~isnan(pQslopes) & ~isinf(pQslopes);
pQslopesM = mean(pQslopes(inds))
pQslopesV = std(pQslopes(inds))
pQslopesVM = pQslopesV/abs(pQslopesM)

figure
subplot(1,3,1)
plot(xE,pE,'k.',xA,pA,'ro',xH,pH,'b*',xQ,pQ,'g+');%xL,pL,'m-s',
xlabel('x')
ylabel('p(x)')
legend({'empirical','parametric','equidistant hist','equiprob hist'})
axis tight;
grid on;
hold off;
subplot(1,3,2)
plot(xElog,pElog,'k.',xAlog,pAlog,'ro',xHlog,pHlog,'b*',xQlog,pQlog,'g+'); %xLlog,pLlog,'m-s',
title(['alpha=',num2str(paramsA(1)),', beta=',num2str(paramsA(2)),', gamma=',num2str(paramsA(3)),', delta=',num2str(paramsA(4))])
xlabel('ln(x)')
ylabel('ln(p(x))')
grid on;
axis equal
hold off;
subplot(1,3,3)
plot(xE,cA,'ro',xE,cE(2:end),'k.'); %xLlog,pLlog,'m-s',
xlabel('x')
ylabel('cdf(x)')
legend({'empirical','parametric'})
axis tight;
grid on;
hold off;
