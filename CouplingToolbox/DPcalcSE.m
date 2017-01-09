function SE = DPcalcSE(Dphi,Nbins,grid)


% Tass P, Rosenblum MG, Weule J, Kurths J, Pikovsky AS, Volkmann J, Schnitzler A, Freund H. 1998. 
% n:m Phase Locking from Noisy Data: Application to Magnetoencephalography. 
% Physical Review Letters. 3291-3294.


% %Wrap phases in the interval (-pi pi)
% Dphi=mod(Dphi,2*pi);
% Dphi(Dphi>pi) = Dphi(Dphi>pi)-2*pi;


%Calculate probability distribution (histogram)
p=hist(Dphi,grid)/length(Dphi);

%Calculate Shannon Entropy
SE = sum(-p(p>0).*log(p(p>0)));

%Normalize with maximum value
norm = log(Nbins);
SE = (norm-SE)/norm;

end