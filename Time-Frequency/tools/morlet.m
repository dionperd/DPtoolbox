function [cwt,scales] = dp_morlet(x,Fs,Freq,w0)

% MORLET CONTINOUS WAVELET TRANSFORM
% [cwt,scales,coi] = morlet(x,Fs,Freq,w0);

  
% Pad signal with zeros to avoid wrap-around
n = length(x);
j = ceil(log2(n));
x2 = zeros(2^j,1);
x2(1:n) = x;

% Transform to Fourier space for fast convolution
N = length(x2);
Xk = fft(x2);

% Make an angular frequency array
k = 1:N;
index1 = find(k <= N/2);
index2 = find(k > N/2);
w = k;
w(index1) = 2*pi*k(index1)*Fs/N;
w(index2) = - 2*pi*k(index2)*Fs/N;

% Automatic scale selection (e.g. for Wavelet denoising)
if (isempty(Freq) == 1)
   Nyquist = Fs/2;
   conversion = (4*pi)/(w0 + sqrt(2 + w0^2));
   s0 = 1./(Nyquist*conversion); % smallest scale = Nyquist
   dj = .25;
   J = log2(N/(Fs*s0))/dj;
   j = 0:J;
   scales = s0 * 2.^(j*dj);
   
% Predefined frequency grid (e.g. for TF analysis)
else
   conversion = (4*pi)/(w0 + sqrt(2 + w0^2));
   scales = 1./(Freq*conversion);
end

scales = scales';

% Compute wavelet transform
m = length(scales);
cwt = ones(m,N);
for trials = 1:m
   s = scales(trials);
   mother = (pi^-0.25).*(w>0).*exp(-(s*w - w0).^2./2);
   daughter = sqrt(2*pi*s*Fs)*mother;
   Wn = ifft(Xk.*daughter');
   cwt(trials,:) = Wn';
end

% Ignore paddes zeros
cwt = cwt(:,1:n); 