function [TF, param] = DPtfviaFFT(x,f,Fs,transform,param,plotting) 

%This function calculates time-frequency transform of several time-series,
%using either the Gabor transform (fixed time and frequency resolution)
%with parameter gamma, or the Gabor (or Morlet) wavelet transform with
%parameter c. The normalization is such that you recover the exact
%amplitude for each time and frequency point.

%This code convolutes by multiplying in the frequency domain (fast).

%
%Inputs:
% -x: the (real) time series stacked in whatever kind of matrix (single or double)
%     where time has to be in the FIRST dimension.
% -f: vector of frequencies, real numbers, where the transformations will
%     be centered
% -Fs: the sampling frequency, real scalar
% -transform: either a string choosing between the up to now available
%            'Gabor' and 'GaborWavelet' (Morlet) transforms, or a set of
%            time-frequency "atoms" (a complex number matrix with dimensions [Nt Nf]) that of some other transform that can be
%            convoluted with x in order to transform it
% -param: optional parameter (gamma for Gabor transformation and c for the
%         Morlet one), it defaults according to the respective comments
%         below
% -plotting: optional flag that if 'true' or positive integer and if the
%            first signal's time-frequency amplitude will be plotted
%            it defaults to 0.
%
%Output:
%  -TF: the complex time-frequency transformed data, a matrix of complex
%       numbers where the first dimension is time, the second frequency,
%       and all the rest follow the 2:end dimensions of x
% -param: the parameter used for the transform (gamma for Gabor transformation and c for the
%         Morlet one)


if nargin<6
    plotting=0;
end

%Supported methods for the moment:
methods = {'gabor','gaborWavelet'};

%Basic relationships:
%Wavelet duration = 2*sigma
%Wavelet length: WL = 2pi*sigma
%Spectral bandwidth = 2*sigmaF
%Morlet parameter c = number of cycles
%c=f/sigmaF = (2pi*sigma)*f = WL*f
%sigma = 1/(2*pi*sigmaF)
%WL = 2*pi*sigma
%where sigma is std in time, sigmaF is std in frequency, c is the cycles
%number and WL is the wavelength.

%Vectorize inputs
xSIZE=size(x);
Nt=xSIZE(1); %the total number of time points
tN=1:Nt;
f=f(:);
f=f.';
Nf=length(f); %the total number of frequency bins

%Reshape x to become a 2D matrix
Nts =prod(xSIZE(2:end)); %total number of time series
x = reshape(x,[Nt, Nts]);

% Ts=1/Fs; %the sampling time
% T = t(end)-t(1); %time length/duration
NFFT = 2^nextpow2(Nt);
f1    = 0:floor(NFFT/2);
f2    = floor(-(NFFT)/2)+1:1:-1;

fall  = [f1 f2]*Fs/NFFT;
fall=fall.';
%Nfall = length(fall);

TwoPi = 2*pi; %a useful constant

%Choose method and construct the Atom of the transform to be convolved with
%each time signal
if ischar(transform) && any(strcmpi(transform,methods))
    
    if strcmpi(transform,'gabor')
        
        if nargin>4
            gamma=param;
        else
            gamma=2*sqrt(pi); %default is the value chosen by Gruber
            param=gamma;
        end
        
        %Gabor function is a Gaussian function with standard deviation:
        sigma = sqrt(1/(2*pi))/gamma;
        %that when gamma=1 allows signals of 
        %at least 2 secs to pass with a weight of at least 4.321% of maximum value
        %at least 1 secs to pass with a weight of at least 45.59% of
        %maximum value,
        %or in other words of a wave length of 
        %WL = 2*pi*sigma = sqrt(2*pi))/gamma = 2.5066/gamma secs
        %which means c ~= 2.5 cycles for the frequency of 1 Hz
        %and a time resolution of 2*sigma = 0.3183 secs
        %and a fequency resolution of 2*sigmaF = 0.3183 Hz as well
        %for all frequencies!
        
        %Correspondance with Gruber's software: 
        %Gaborfenster     = sqrt(1/(2*pi)) * exp(-(Frequenzen - om).^2 ./(4*MD.Gamma));
        
        %In frequency space the standard deviation becomes:
        sigmaF = 1/(2*pi*sigma);
        
        
        FALL = repmat(fall,[1,Nf]);
        F = repmat(f,[NFFT,1]);
        AtomF = exp( -(FALL-F).^2 / (2*sigmaF^2) );
        
       %Correspondance with Gruber's software: 
       %Gaborfenster     = sqrt(1/(2*pi)) * exp(-(Frequenzen - om).^2./(4*MD.Gamma));
       %i.e. 2*sigmaF^2 = 4*MD.Gamma => sigmaF = sqrt(2*MD.Gamma) 
       %=> sigmaF = 1/(2*pi*sqrt(2*MD.Gamma)),
       %so that when MD.Gamma = 1, we have 
       %sigmaF = sqrt(2) = 1.4142 Hz 
       %and sigma = sqrt(2)/(4*pi) = 0.1125 sec, 
       %& a wave length of WL = 2*pi*sqrt(2)/(4*pi) = 1/sqrt(2) = 0.7071 sec
       %or c ~= 0.7071 cycles at the frequency of 1 Hz,
       %and a time resolution of 2*sigma = 0.2251 secs
       %and a fequency resolution of 2*sigmaF = 2.8284 Hz 
       %for all frequencies!
       
       %so, to get the same hear we need gamma = 2*sqrt(pi)
        
    elseif strcmpi(transform,'gaborWavelet')
        
        
        
        %Default morlet parameter: c=7 cycles per wavelet
        %which leads to a standard deviation of sigma = 7/(2*pi*f)
        %so that for f=1Hz 
        %we have a wave length of WL = 2*pi*sigma =7 secs 
        %a time resolution of 2*sigma = 2.2282 secs
        %and a frequency resolution of 2*sigmaF = 0.2857 Hz
        %and for f=50 Hz, we have WL = 0.8796
        %a time resolution of 0.0446 sec
        %and a frequency resolution of 14.2857 Hz
        if nargin>4
            c=param;
        else
            c=7;
            param=c;
        end
        sigma = @(f)max(c./(TwoPi*f),eps);
        %In frequency space the standard deviation becomes:
        sigmaF = 1./(2*pi*sigma(f));
        
        FALL = repmat(fall,[1,Nf]);
        F = repmat(f,[NFFT,1]);
        SIGMAF = repmat(sigmaF,[NFFT,1]);
        
        AtomF = exp(-(FALL-F).^2./(2*SIGMAF.^2));

        %Lachaux et al 1999 use sigma = 7/f => sigmaF = f/(14*pi) => c=14*pi
        %so that for f=1Hz we have a wavelength of WL = 2*pi*sigma = 43.9823 secs
        %a time resolution of 2*sigma = 14 secs
        %and a frequency resolution of 0.0455 Hz
        %and for f=50Hz, WL = 2*pi*sigma =  0.8796 secs  
        %time resolution 0.28 sec
        %and frequency resolution 2.2736 Hz
    end
    
    
elseif isnumeric(transform) && all(size(transform)==[Nt Nf])
    
    AtomF = transform;
    
else
    error('transform is neither a string for one of the supported methods nor a matrix of size [Ntimes Nfreqs]')
end


%Execute the transform
if isa(x, 'single')
    TF = nan(Nt,Nf,Nts,'single');
else
    TF = nan(Nt,Nf,Nts);
end
for iF = 1:Nf;
    temp=2*ifft( repmat( AtomF(:,iF), [1,Nts] ).*fft(x,NFFT,1) );
    TF(tN,iF,:)= temp(tN,:);
end

% %Older code that wouldn't allowed arbitrary dimensions for x
% Nsiz = length(xSIZE);
% Ns1 = xSIZE(2);
% switch Nsiz
%     
%     case 2
%         
%         if Ns1>1
%             if isa(x, 'single')
%                 TF = nan(Nt,Nf,Ns1,'single');
%             else
%                 TF = nan(Nt,Nf,Ns1);
%             end
%             for iF = 1:Nf;
%                 temp=2*ifft( repmat( AtomF(:,iF), [1,Ns1] ).*fft(x,NFFT,1) );
%                 TF(tN,iF,:)= temp(tN,:);
%             end
%                
%         else
%             if isa(x, 'single')
%                 TF=nan(Nt,Nf,'single');
%             else
%                 TF = nan(Nt,Nf);
%             end
%             for iF = 1:Nf;
%                 temp = 2*ifft(  AtomF(:,iF).*fft(x,NFFT) );
%                 TF(tN,iF)= temp(tN);
%             end
%         end
%         
%     case 3
%         Ns2=xSIZE(3);
%         if isa(x, 'single')
%             TF = nan(Nt,Nf,Ns1,Ns2,'single');
%         else
%             TF = nan(Nt,Nf,Ns1,Ns2);
%         end
%         for iF = 1:Nf;
%             temp = 2*ifft( repmat( AtomF(:,iF), [1,Ns1,Ns2] ).*fft(x,NFFT,1) );
%             TF(tN,iF,:,:)= temp(tN,:,:);
%         end
%         
%     otherwise
%         error('x has a size bigger than the supported ones')
% end


%Plotting if required...
%(it will plot the first signal only, just for demonstration...)
if (plotting)
    t=[1:Nt]*1/Fs;
    %Time-frequency analysis
    figure;
    subplot(2,1,1)
    imagesc(t,f,abs(TF(:,:,1)).');
    set(gca,'Ydir', 'normal')
    xlabel('Time (sec)');
    ylabel('Frequency (sec)');
    hold off;
    subplot(2,1,2)
    imagesc(t,f,angle(TF(:,:,1)).');
    set(gca,'Ydir', 'normal')
    xlabel('Time (sec)');
    ylabel('Frequency (sec)');
    hold off;
    % close(h);
    % clear h;
end


%Reshape TF into the original size of x:
TF = squeeze(reshape(TF,[Nt,Nf,xSIZE(2:end)]));

