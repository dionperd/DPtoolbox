function PLI=DPcalcPLI(Dphi)


% Stam CJ, Nolte G, Daffertshofer A. 2007. 
% Phase lag index: assessment of functional connectivity from multi channel EEG and MEG with diminished bias from common sources. 
% Human brain mapping. 28:1178?1193.

% %Wrap phases in the interval (-pi pi)
% Dphi=mod(Dphi,2*pi);
% Dphi(Dphi>pi) = Dphi(Dphi>pi)-2*pi;


%Calculation of phase lag index
PLI = mean(sign(Dphi));