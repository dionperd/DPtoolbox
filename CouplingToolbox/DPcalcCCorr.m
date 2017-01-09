function CCR=DPcalcCCorr(phi)% 

% 
% Burgess AP. 2013. 
% On the interpretation of synchronization in EEG hyperscanning studies: a cautionary note. 
% Front Hum Neurosci. 7:881.
%There is a mistake in equation 6 of this paper.
%The correct equation is 
%CCorr = sum(sin(phi1-phi1_mean)*sin(phi2-phi2_mean)) / 
%        sqrt(sum(sin(phi1-phi1_mean).^2)*sum(sin(phi2-phi2_mean).^2))

% %Wrap phases in the interval (-pi pi)
% phi=mod(phi,2*pi);
% phi(phi>pi,:) = phi(phi>pi,:)-2*pi;


%deviations from the mean angle for each point
phiM = mean(phi); %the mean phases for each signal
x1 = sin(phi(:,1)-phiM(1));
x2 = sin(phi(:,2)-phiM(2));

nom = sum(x1.*x2);
den = sqrt(sum(x1.^2).*sum(x2.^2));

%Calculation of circular correlation coefficient
CCR=nom/den;


end