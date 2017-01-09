function [PLV]=DPcalcPLV(Dphi)% PP SG


% Lachaux J-P, Rodriguez E, Le van Quyen M, Lutz A, Martinerie J, Varela FJ. 2000. 
% Studying single-trials of phase synchronous activity in the brain. 
% International Journal of Bifurcation and Chaos. 10:2429?2439.
% 
% Lachaux J-P, Rodriguez E, Martinerie J, Varela FJ. 1999. 
% Measuring phase synchrony in brain signals. 
% Human brain mapping. 8:194?208.


% %Wrap phases in the interval (-pi pi)
% Dphi=mod(Dphi,2*pi);
% Dphi(Dphi>pi) = Dphi(Dphi>pi)-2*pi;


%Calculation of phase locking value
PLV=abs( mean( exp( 1i * Dphi  ) ) );



end