function x=DPsurrShufflPhase(x,D,Ntr,NFFT)

%This function shuffles randomly the phases of the columns of x in Fourier
%space in order to create surrogate time series

%IMPORTANT: in order that the inverse Fourier transform produces real
%numbers, FFTx has to be conjugate symmetric,
%i.e., FFTx(i,:) = conj(FFTx(mod(N-i+1,N)+1,:))

%which means that N (number of time points) has to be odd number!!!!

%Get the initial size of x:
sizX0 = size(x);
reshapeFlag = length(sizX0)>2;
N=sizX0(1);
if nargin<4
    NFFT=N;
end
if reshapeFlag
    %Reshape so that only the first (time) dimension is retained and all
    %others are stacked to the second one
    x=reshape(x,[N,prod(sizX0(2:end))]);
end

%Calculate amplitude and phase of fft
X=fft(x);
absX = abs(X);
anglX = angle(X);
clear X;

%Get the middle point
N2=floor(N/2);
ii=[1:N2]+1; 

%Permute the first half of the phase time series, excluding the first
%element (Frequency = 0, constant term)

iiRand = randperm(N2)+1;
anglX(ii,:) =  anglX(iiRand,:);

%Create the new complex signal in Fourier space for this half with
%amplitude 1
X(ii,:) = exp(1i*anglX(ii,:));

%Calculate the conjugate (second half) part
ii=ii(end)+1:N;
X(ii,:) = conj(X(mod(N-ii+1,N)+1,:));

%Multiply with the initial amplitude:
X=X.*absX;

%Calculate the new -phase scarambled- signal in time domain with inverse
%Fourier transform
x=ifft(X,NFFT);
if ~isreal(x)
    disp('WARNING: Phase scrambled signal is not real! Minimum and maximum imag parts are:')
    imagX = imag(x);
    min(imagX(:))
    max(imagX(:))
    x=real(x);
    disp('Consider discarding one data point so that they are of an odd number')
end
    
if reshapeFlag
    %Reshape to the original size:
    x=reshape(x,sizX0);
end
