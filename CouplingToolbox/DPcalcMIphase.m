function MI = DPcalcMIphase(phi,Nbins,grid)

%Similar to:
% Tass P, Rosenblum MG, Weule J, Kurths J, Pikovsky AS, Volkmann J, Schnitzler A, Freund H. 1998. 
% n:m Phase Locking from Noisy Data: Application to Magnetoencephalography. 
% Physical Review Letters. 3291-3294.

N=size(phi,1);


%For each signal... 
for iS=1:2;
        
    %...calculate probability distribution (histogram)...
    p=hist(phi(:,iS),grid)/N;
    
    %...calculate Shannon entropy
    H(iS) =  sum(-p(p>0).*log(p(p>0))) / log(Nbins);
    
end

%Calculate the joint distribution...
p=hist3(phi,{grid grid})/N;
p=p(:); %(vectorize)
%...and the joint Shannon entropy
H12 = sum(-p(p>0).*log(p(p>0))) / log(Nbins);

%Calculate mutual information:
MI = sum(H) - H12;


end