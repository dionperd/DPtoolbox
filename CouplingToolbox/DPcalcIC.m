function [PCI NCI ACI ICI12 ICI21] = DPcalcIC(Dphi,T,Dphi0)

% Müller V, Lindenberger U. 2011. 
% Cardiac and respiratory patterns synchronize between persons during choir singing. 
% PloS one. 6:e24893.
% 
% Kitzbichler MG, Smith ML, Christensen SR, Bullmore E. 2009. 
% Broadband criticality of human brain network synchronization. 
% PLoS computational biology. 5:e1000314.


% %Inputs:
% %-Dphi: phase differences, Dphi~[-pi,pi], real numbers in a vector format
% %-Dphi0: phase difference threshold for synchronization Dphi~(0,pi/2)
% %-T: period in time points, positive, integer, 0<T<length(Dphi)
% 
% %Outputs:
% %....The measures PCI, NCI, ACI, ICI, real numbers ~[0,1]
% 
% 
% %Data validation
% 
% funcName = 'DP_PCI_NCI_ACI_ICI';
% 
% varName='Dphi';
% testDPHI = { @(Dphi)isnumeric(Dphi),...
%              @(Dphi)isreal(Dphi),...
%              @(Dphi)min(Dphi)>=-pi,...
%              @(Dphi)max(Dphi)<=pi,...
%              @(Dphi)isvector(Dphi),...
%              @(Dphi)size(Dphi,2)==1 };
% param={{},{},{},{},{},{}};
% mode=['e','e','e','e','e','w'];
% execfun={{},{},{},{},{},@(Dphi)Dphi.'};
% default=nan;
% [~, Dphi] = DPvalidateData(Dphi,testDPHI,param,mode,execfun,default,varName,funcName);
% 
Nphi = length(Dphi); %The total number of points
% 
% 
% varName='Dphi0';
% testDPHI0 = { @(Dphi0)isnumeric(Dphi0),...
%               @(Dphi0)isreal(Dphi0),...
%               @(Dphi0)Dphi0>=0,...
%               @(Dphi0)Dphi0<=pi/2,...
%               @(Dphi0)isscalar(Dphi0) };
% param={{},{},{},{},{}};
% mode=['e','e','e','e','e'];
% execfun={{},{},{},{},{}};
% default=pi/4;
% [~, Dphi0] = DPvalidateData(Dphi0,testDPHI0,param,mode,execfun,default,varName,funcName);
% 
% 
% varName='T';
% testT = { @(T)isnumeric(T),...
%           @(T)isreal(T),...
%           @(T)T>0,...
%           @(T,Nphi)T<Nphi,...
%           @(T)round(T)==T,...
%           @(T)isscalar(T) };
% param={{},{},{},{Nphi},{},{}};
% mode=['e','e','e','e','e','e'];
% execfun={{},{},{},{},{},{}};
% default=nan;
% [~, T] = DPvalidateData(T,testT,param,mode,execfun,default,varName,funcName);


% %Wrap phases in the interval (-pi pi)
% Dphi=mod(Dphi,2*pi);
% Dphi(Dphi>pi) = Dphi(Dphi>pi)-2*pi;

%If the time window is smaller than a period...
if Nphi<T
    %...adjust...
    T=Nphi;
    %...and print a warning message:
    cprintf('Magenta',['WARNING: Time window is smaller than one period. Clearing period is adjusted to NwinLen=%d time points.\n'],T);
end

%Initialize the meters for...
Na=0; %...absolute sync points
Np=0; %...positive sync points
Nn=0; %...negative sync points

%Initialize the temporary meters for...
tempNa=0; %...absolute sync points
tempNp=0; %...positive sync points
tempNn=0; %...negative sync points

%Add a phase difference greater than the threshold at the end, so that we
%stop with a non sync point
Dphi = [Dphi;Dphi0+pi/10]; 
         
%Main loop
%For every point in the sequence...
for ii=1:Nphi+1;
    
    %...check if it is an absolute sync point
    if abs(Dphi(ii))<=Dphi0
        
        %...if yes, increase tempNa...
        tempNa=tempNa+1;
        
        %...and check if it is positive
        if Dphi(ii)>=0
            %...if yes, increase tempNp...
            tempNp=tempNp+1;
        else %...if not, increase tempNn...
            tempNn=tempNn+1;
        end
    else %...if it is not a sync point
        %...check if we have equal or more consecutive sync points than 1 period
        if tempNa>=T
            %...if yes, count the temporary sync points as sync points
            Na=Na+tempNa;
            Np=Np+tempNp;
            Nn=Nn+tempNn;
        end
        
        %...in any case, reset all temporary meters to 0
        tempNp=0;
        tempNn=0;
        tempNa=0;       
    end     
end

%Finally, calculate measures according to their definitions.
PCI = Np/Nphi;
NCI = Nn/Nphi;
ACI = Na/Nphi;
if ACI~=0
    ICI12 = ( (PCI+ACI)/(2*ACI) )*sqrt(PCI);
    ICI21 = ( (NCI+ACI)/(2*ACI) )*sqrt(NCI);
else
    ICI12=0;
    ICI21=0;
end

