function [spectrum,freqoi,timeoi] = ft_jakob_gabor(dat, time, varargin)
% [nchan,nPoints] = size(dat);
%  gamma = 1;
%  TapBreite = .15;
%  
% srate = 250;
% deltaf  = srate/nPoints;          %Frequenzauflösung
% 
% T       =   nPoints / srate; % Zeit in Sekunden
% f1         = 0:floor(nPoints/2);
% f2         = floor(-(nPoints)/2)+1:1:-1;
% f          = [f1 f2];
% Frequenzen = f * deltaf;
% 
% 
% frequency_legend = Frequenzen >= 1 &  Frequenzen <=40;
% frequency_legend = Frequenzen(frequency_legend);
% 
% AnzFreq    = length(frequency_legend);%ceil((frequency_interval(2)-frequency_interval(1)) / deltaf);
% TapWin  = tukeywin(nPoints,TapBreite)';
% spectrum = nan(nchan,AnzFreq,nPoints);
% for iChan = 1:nchan
%     data = TapWin.*dat(iChan,:);
%     
%     Y        = (1/(sqrt(T) * srate)) * fft(data);           %first transform into frequency s
%     Transformed = nan(AnzFreq,nPoints);
%     for iFreq = 1 : AnzFreq;
%         om = frequency_legend(iFreq);
%         Gaborfenster     = sqrt(1/(2*pi)) * exp(-(Frequenzen - om).^2 ./ (4*gamma));
%         sumGaborfenster=sum(Gaborfenster);
%         Gaborfenster = Gaborfenster/sumGaborfenster; %make the area under the curve = 1
%         Transformed(iFreq,:) = ifft( Y .* Gaborfenster);
%     end;
%     spectrum(iChan,:,:) = Transformed;
% end
% disp('')
% 
% freqoi = frequency_legend;
% timeoi = time;
% return

freqoi      = ft_getopt(varargin, 'freqoi', 'all');
timeoi      = ft_getopt(varargin, 'timeoi', 'all');
pad       = ft_getopt(varargin, 'pad');
polyorder = ft_getopt(varargin, 'polyorder', 0);

gamma       = ft_getopt(varargin, 'gamma', 1);
TapBreite   = ft_getopt(varargin, 'TapBreite', .15);


%% this comes from ft_specest_wavelet
% Set n's
[nchan,ndatsample] = size(dat);

% Remove polynomial fit from the data -> default is demeaning
% if polyorder >= 0
%   dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
% end

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad = zeros(1,round((pad - dattime) * fsample));
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data

% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all') % if input was 'all'
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
% check for freqoi = 0 and remove it, there is no wavelet for freqoi = 0
if freqoi(1)==0
  freqoi(1)  = [];
  freqboi(1) = [];
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);

% Set timeboi and timeoi
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
  timeoi   = round(timeoi .* fsample) ./ fsample;
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

%% now the gabor transform
disp('')

TapWin  = tukeywin(ndatsample,TapBreite)';
dat = dat .* repmat(TapWin,[nchan,1]);

% from Gruber's EEGnew
deltaf      = fsample/ndatsample; 
T           = ndatsample / fsample; % Zeit in Sekunden
f1         = 0:floor(ndatsample/2);
f2         = floor(-(ndatsample)/2)+1:1:-1;
f          = [f1 f2];
Frequenzen = f * deltaf;

Y        = (1/(sqrt(T) * fsample)) * fft(dat,[],2);           %first transform into frequency space

% tictoc = tic;
% spectrum = nan(nchan,nfreqoi,ndatsample);
% for iFreq = 1 : nfreqoi;
%     om = freqoi(iFreq);
%     Gaborfenster     = sqrt(1/(2*pi)) * exp(-(Frequenzen - om).^2 ./ (4*gamma));
%     Gaborfenster = Gaborfenster/sum(Gaborfenster); %make the area under the curve = 1
%     spectrum(:,iFreq,:) = ifft( Y .* repmat(Gaborfenster,[nchan, 1]),[],2);
% end
% tictocA = toc(tictoc);

% get rid of the for-loop to improve performace
Gaborfenster = nan(nfreqoi,ndatsample);
for iFreq = 1 : nfreqoi
    om = freqoi(iFreq);
    Gaborfenster(iFreq,:)     = sqrt(1/(2*pi)) * exp(-(Frequenzen - om).^2 ./ (4*gamma));
end

Gaborfenster = Gaborfenster./repmat(sum(Gaborfenster,2),[1 ndatsample]); %make the area under the curve = 1
Gaborfenster = repmat(permute(Gaborfenster,[3,1,2]),[nchan, 1,1]);
Y = repmat(permute(Y,[1,3,2]),[1,nfreqoi,1]);
spectrum = ifft( Y .* Gaborfenster,[],3);

%only take the timepoints of interest
spectrum = spectrum(:,:,timeboi);
