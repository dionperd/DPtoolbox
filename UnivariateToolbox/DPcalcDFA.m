function [H, logF ] = DPcalcDFA(DATA,scales,logScales,order,N,Nsc,slopeInds) %DP: plotting

% DFA performs the Detrended Fluctuations Analysis of a given time series as a
% function of scales and delivers it as output1 of this function.
%H: Hurst exponent
%F: log of fitting at each scale

% DATA = timeseries
% win_length = desired window length of a given box, needs to be specified in fucntion  
% order = order of polynomial fit, for DFA to be set as 1
%
% 10/04/2011 - Original version: Viktor Jirsa
%              this version used the functions LinearFit and LinearLSQ
% 03/05/2011 - VJ: simplified version of the algorithm using the MATLAB
%              functions polyfit and polyval (found on the internet, tested
%              by VJ)

%DP modified it to loop over calculations for different scales and finally 
%calculate the slope of the log(scale)-log(F) line, as well as plot, 22-04-2014

%N=length(DATA);        % total length of time series

%DP added the loop for the different scales
%Nsc = length(scales);


F=zeros(Nsc,1);
for iSc=1:Nsc;
    win_length = scales(iSc);
    
    if win_length>(order + 1)
    
        WL=[1:win_length].';
         
        n=floor(N/win_length); % number of boxes on a given scale
        N1=n*win_length;
        y=zeros(N1,1); Yn=zeros(N1,1); fitcoef=zeros(n,order+1);
        
        mean1=mean(DATA(1:N1));
        for i=1:N1
            y(i)=sum(DATA(1:i)-mean1);
        end

        for j=1:n
            fitcoef(j,:)=polyfit(WL,y(((j-1)*win_length+1):j*win_length),order);
            Yn(((j-1)*win_length+1):j*win_length)=polyval(fitcoef(j,:),WL);
        end
        
        output1=sqrt(sum((y-Yn).^2)/N1);
        
        F(iSc) =output1;
        
    else
        F(iSc) = nan;
    end
end

%logScales = log(scales);
logF = log(F);

% H=zeros(Nsc,1);
% inF = ~isinf(logF);
% % p=polyfit( logScales(inF), logF(inF), 1);
% % H = p(1);
% H(inF) = DPcalcSlope(logScales(inF), logF(inF), sum(inF), Nslope);

%H=zeros(Nsc,1);
% inF = ~isinf(logF);
% p=polyfit( logScales(inF), logF(inF), 1);
% H = p(1);
% H(inF) = DPcalcSlope(logScales(inF), logF(inF), sum(inF), Nslope);

if ~isempty(slopeInds)
%     logScales=logScales(slopeInds);
%     logFtemp=logF(slopeInds);
%     inF = ~isinf(logFtemp) & ~isnan(logFtemp);
%     p=polyfit( logScales(inF), logFtemp(inF), 1);
%     H=p(1);
    H = DPcalcSlope(logScales,logF,slopeInds);
else
    H=0;
end
logF = logF.';
% if plotting
%     plot(logScales,logF)
% end


