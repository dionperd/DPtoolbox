function powORampl = DPcalcPowerORampl(x,fs,band,cutTails,mode)


% Cohen MX. 2008. Assessing transient cross-frequency coupling in EEG data. 
% Journal of Neuroscience Methods. 168:494?499.


%x: (time x trials)

%Calculate the amplitude or power of signal x
[~,powORampl,~]= DPco_hilbproto(x,cutTails);
if strcmpi(mode,'pow')
    powORampl=powORampl.^2;
end
%Detrend it by subtracking the mean
N=size(powORampl,1);
powORampl = powORampl - repmat(mean(powORampl,1),N,1);

%Filter the power with a bandpass filter
%[smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff,epochframes,filtorder,revfilt,firtype)
%                                    epochframes,filtorder,revfilt,firtype);
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high-edge frequency in pass band (Hz) {0 -> highpass}
%   epochframes = frames per epoch (filter each epoch separately {def/0: data is 1 epoch}
%   filtorder   = length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {0}
%   firtype     = 'firls'|'fir1' {'firls'}
%
% Outputs:
%    smoothdata = smoothed data
%    filtwts    = filter coefficients [smoothdata <-
%    filtfilt(filtwts,1,data)]
powORampl=DPeegfilt(powORampl.',fs,band(1),band(2));
powORampl=powORampl.';
