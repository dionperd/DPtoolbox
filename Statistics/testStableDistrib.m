%testStableDistrib


%%Create distribution

%parameters
alpha = 0.5;
beta = 1;
gamma = 1;
delta = 1;

%support
x = [-10:0.1:10].';
indP = x>0;
xlog = log(x(indP));

%number of samples
N = 1000;

%alpha stable:
tic
xA = stblrnd(alpha,beta,gamma,delta,[N,1]);
pA = stblpdf(x,alpha,beta,gamma,delta);
pAlog = log(pA(indP));
toc

%LSD:
tic
xL = levystblrnd(alpha, beta, gamma, delta,[N,1]);
pL=levystblpdf(x, alpha, beta, gamma, delta);
pLlog = log(pL(indP));
toc

% %Plot distributions
% figure;
% subplot(1,2,1)
% plot(x,pA,'r');
% hold on;
% plot(x,pL,'m');
% grid on;
% hold off;
% subplot(1,2,2)
% plot(xlog,pAlog,'r');
% hold on;
% plot(xlog,pLlog,'m');
% grid on;
% hold off;


%Estimate distributions

%parametric estimation

%alpha stable:
%fit
tic
paramsAA = stblfit(xA,'ecf',statset('Display','iter'));
paramsLA = stblfit(xL,'ecf',statset('Display','iter'));
toc
paramsAA
paramsLA
%calculate
pAA = stblpdf(x,paramsAA(1),paramsAA(2),paramsAA(3),paramsAA(4));
pLA = stblpdf(x,paramsLA(1),paramsLA(2),paramsLA(3),paramsLA(4));
pAAlog = log(pAA(indP));
pLAlog = log(pLA(indP));

%LSD:
%fit
tic
paramsAL = levystblfit(xA);
paramsLL = levystblfit(xL);
toc
paramsAL
paramsLL
%calculate
pAL = levystblpdf(x,paramsAL(1),paramsAL(2),paramsAL(3),paramsAL(4));
pLL = levystblpdf(x,paramsLL(1),paramsLL(2),paramsLL(3),paramsLL(4));
pALlog = log(pAL(indP));
pLLlog = log(pLL(indP));


%non-parametric estimation
Nbins = sqrt(N);

%equidistang histogram calculation
tic
[pAh, xAh, xiAh] = DPcalcPDF_1D(xA,'equidist',Nbins,0);
[pLh, xLh, xiLh] = DPcalcPDF_1D(xL,'equidist',Nbins,0);
toc
temp = xAh>0;
xAhlog = log(xAh(temp));
pAhlog = log(pAh(temp));
temp = xLh>0;
xLhlog = log(xLh(temp));
pLhlog = log(pLh(temp));


%equidistang histogram calculation
tic
[pAq, xAq, xiAq] = DPcalcPDF_1D(xA,'isohisto',Nbins,0);
[pLq, xLq, xiLq] = DPcalcPDF_1D(xL,'isohisto',Nbins,0);
toc
temp = xAq>0;
xAqlog = log(xAq(temp));
pAqlog = log(pAq(temp));
temp = xLq>0;
xLqlog = log(xLq(temp));
pLqlog = log(pLq(temp));


figure
subplot(2,2,1)
plot(x,pA,'r--',x,pAA,'r-o',x,pAL,'m-o',xAh,pAh,'b',xAq,pAq,'g');
grid on;
hold off;
subplot(2,2,2)
plot(x,pL,'m',x,pLL,'m-o',x,pLA,'r-o',xLh,pLh,'b',xLq,pLq,'g');
grid on;
hold off;
subplot(2,2,3)
plot(xlog,pAlog,'r--',xlog,pAAlog,'r-o',xlog,pALlog,'m-o',xAhlog,pAhlog,'b',xAqlog,pAqlog,'g');
grid on;
hold off;
subplot(2,2,4)
plot(xlog,pLlog,'m',xlog,pLLlog,'m-o',xlog,pLAlog,'r-o',xLhlog,pLhlog,'b',xLqlog,pLqlog,'g');
grid on;
hold off;
