function [FFT, f, t] = DP_fft(x,Fs)

SIZ=size(x);
if SIZ(2)>SIZ(1)
    x=x';
end

L=length(x);
Ts=1/Fs;
t = (0:L-1)*Ts;
NFFT = 2^nextpow2(L);
FFT = fft(x,NFFT,1)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
plot(f,2*abs(FFT(1:NFFT/2+1,:)));grid on;
xlabel('Hz')
ylabel('Amplitude')
hold off
