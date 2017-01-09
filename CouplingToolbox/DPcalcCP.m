function CP = DPcalcCP(phi,Nbins,grid)

% Tass P, Rosenblum MG, Weule J, Kurths J, Pikovsky AS, Volkmann J, Schnitzler A, Freund H. 1998. 
% n:m Phase Locking from Noisy Data: Application to Magnetoencephalography. 
% Physical Review Letters. 3291?3294.

%Loop for the measure calculation
CP=0;

%For each bin...
for iB=1:Nbins;
    
    
    %Find the indexes of this phi phases that lie in this bin
    thisPhi1BinInds =  (phi(:,1)>=grid(iB))&(phi(:,1)<grid(iB+1));
    
    %Get the subset of the other phi phases that correspond to the above phi1 phases
    %(devide my n (for phi2) or m (for phi1) respectively so that the
    %range of the phase distribution is identical)
    thisPhi2Bin = phi(thisPhi1BinInds,2);
    
    %Add to the measure the value of this bin
    if ~isempty(thisPhi2Bin)
        CP = CP + abs( mean( exp( 1i * thisPhi2Bin ) ) );
    end
end


%Take the average of all bins
CP=CP/Nbins;



end