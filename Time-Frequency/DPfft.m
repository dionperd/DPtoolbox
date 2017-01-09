function [FFT, f] = DPfft(x,Fs,NFFT)

N=size(x,1);

if nargin<3
    NFFT = 2^nextpow2(N);
end

f = Fs/2*linspace(0,1,NFFT/2+1);

FFT = fft(x,NFFT);

subplot(2,1,1)
hold on;
plot(f,2*abs(FFT(1:NFFT/2+1,:)));grid on;
xlabel('Hz')
ylabel('Amplitude')
hold off

subplot(2,1,2)
hold on;
plot(log(f),log(2*abs(FFT(1:NFFT/2+1,:))));grid on;
axis equal
xlabel('log(Hz)')
ylabel('log(Amplitude)')
hold off
