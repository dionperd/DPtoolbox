function PPC=DPcalcPPC(Dphi)

% Vinck M, van Wingerden M, Womelsdorf T, Fries P, Pennartz CM a. 2010. 
% The pairwise phase consistency: a bias-free measure of rhythmic neuronal synchronization. 
% NeuroImage. 51:112?122.


% %Wrap phases in the interval (-pi pi)
% Dphi=mod(Dphi,2*pi);
% Dphi(Dphi>pi) = Dphi(Dphi>pi)-2*pi;

%------------------------------Slow way------------------------------------
% %Find all unique pairs:
% % 
% %1. Create all possible pairs
% [Dphi11 Dphi22]=meshgrid(Dphi); %phase differences
% [i j]=meshgrid(1:length(Dphi)); %indexes
% 
% %2.Keep only pairs where the index of the 1st element 
% %is bigger than the index of the second one,
% %i.e. i>j, where i: index of rows, j: index of columns
% Dphi1=Dphi11(i>j);
% Dphi2=Dphi22(i>j);
% 
% %Calculate now the measure
% PPC= mean( cos(Dphi1-Dphi2) );
%------------------------------Slow way------------------------------------

N=length(Dphi);
PPC=0;
for i=1:N;
    for j=i+1:N;
        PPC = PPC + cos(Dphi(i)-Dphi(j));
    end
end
N=(N*(N-1))/2;
PPC = PPC/N;

end
