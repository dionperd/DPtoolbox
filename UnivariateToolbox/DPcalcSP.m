function [P, DOF, logP, PS] = DPcalcSP(x,logf,winfun,NFFT,Nf,N,slopeInds) %DP: plotting


%This function calculates spectral measures:
%DOF: degrees of freedom, ->0 for sharp peak, ->1 for white noise
%PS: power slope in a logP-logf diagram
%logP: is the log power

% Vaillancourt DE, Newell KM. 2003. 
% Aging and the time and frequency structure of force output variability. 
% J Appl Physiol. 94:903?912.


%These are calculated outside the function for computation speed:
% N=size(x,1);
% NFFT = 2^nextpow2(N);
% Nf=NFFT/2;
% f = fs/2*linspace(0,1,Nf+1).';
% logf=log(f(2:end));



%Calculate fft 
x=diag(winfun(N))*x; %apply the window to all columns of x
FT = fft(x,NFFT);
FT=FT(2:Nf+1,:); %get the positive side of the spectrum, without the DC term
P = FT.*conj(FT); %calculate power
%%pxx=pwelch(x,window,noverlap,nfft,fs)
% P=pwelch(x,winfun(N),0,NFFT);
% P=P(2:end);

%DOF normalized with the number of frequency bins:
DOF = sum(P)^2./sum(P.^2)/Nf;

logP = log(P);

% if Nslope
%     %Calculate logP and fit a line to get the slope:
%     %logf = log(f);
%     
%     %     p=polyfit( logf, logP, 1);
%     %     PS = p(1);
%     PS = DPcalcSlope(logf, logP, Nf, Nslope);
% else
%     PS=zeros(Nf,1);
% end

if ~isempty(slopeInds)
%     logf=logf(slopeInds);
%     logPtemp=logP(slopeInds);
%     inF = ~isinf(logPtemp) & ~isnan(logPtemp);
%     p=polyfit( logf(inF), logPtemp(inF), 1);
%     PS=p(1);
    PS = DPcalcSlope(logf,logP,slopeInds);
else
   PS=0;
end
logP = logP.';

%plot(logf,logP)

