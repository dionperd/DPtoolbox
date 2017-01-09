function [GBT, f, t] = DPgabor(x,fs,t,f,gamma) 

%Vectorize inputs
t=t(:);
f=f(:);
f=f.';
xSIZE=size(x);
Nt = xSIZE(1);       
Nf=length(f);%the total number of frequency bins

Ts=1/fs; %the sampling time

if nargin<5
    gamma=1;
end

%the Gabor window in time as a gaussian window (f(x) = exp(-x^2/(2*c^2)))  
%with c=gamma/(Nt*Ts*sqrt(pi/2))
GaborWin = gausswin(Nt,Nt*Ts*sqrt(pi/2)/gamma); 
GaborWin = 2*GaborWin/sum(GaborWin);%normalize
%Create matrices for faster calculation of the convolutions
GW = repmat(GaborWin,1,Nf);
T = repmat(t,1,Nf);
F = repmat(f,Nt,1);
%The Gabor atom to be convolved with each time signal
GaborAtom = GW.* exp(1i*2*pi*F.*T);

Nsiz = length(xSIZE);
Ns1 = xSIZE(2);
switch Nsiz
    
    case 2
        
        if Ns1>1
            
            GBT = nan(Nt,Nf,Ns1);
            for iS1=1:Ns1;
                GBT(:,:,iS1)=conv2(GaborAtom, x(:,iS1) ,'same');
            end
            
        else
            GBT = nan(Nt,Nf);
            GBT(:,:)=conv2(GaborAtom, x,'same');
        end
        
    case 3
        Ns2=xSIZE(3);
        GBT = nan(Nt,Nf,Ns1,Ns2);
        for iS1=1:Ns1;
            for iS2=1:Ns2;
                GBT(:,:,iS1,iS2)=conv2(GaborAtom, x(:,iS1,iS2) ,'same');
            end
        end
        
end