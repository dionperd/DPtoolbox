function phi = DPpower2phase(x,fs,band,cutTails)


% Cohen MX. 2008. Assessing transient cross-frequency coupling in EEG data. 
% Journal of Neuroscience Methods. 168:494?499.


%x: (time x trials)

powORampl = DPcalcPowerORampl(x,fs,band,[0 0]);

%Calculate the phase of the power in the specific frequency band
[~,phi,~]= DPco_hilbproto(powORampl,cutTails);

end