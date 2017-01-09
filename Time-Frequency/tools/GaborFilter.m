function [Transformed, Spektrum] = GaborFilter(MD, x)

Gamma   = MD.Gamma;
SR      = MD.D_SR;
T       = length(x) / SR;
UBF     = MD.UntereFreq;
OBF     = MD.ObereFreq;
df      = MD.deltaf;
N       = length(x);                     %Anzahl Messwerte
%---------------------- Fouriertransformation des Signals
Y        = (1/(sqrt(T) * SR)) * fft(x);
Spektrum = sqrt(2 * (Y .* conj(Y)));
%---------------------- Fouriertransformation des Signals
AnzFreq    = round((OBF - UBF) / df) + 1;
f1         = 0:floor(N/2);
f2         = floor(-(N)/2)+1:1:-1;
f          = [f1 f2];
Frequenzen = f * df;
om         = ones(1,AnzFreq);
X          = 1:N;
for i = 1 : AnzFreq;
    %om               = i * df; %diskrete Frequenzen, durch erwünschte Auflösung vorgegeben 
    om               = i * df + (UBF); %diskrete Frequenzen, durch erwünschte Auflösung vorgegeben 
    Gaborfenster     = sqrt(1/(2*pi)) * exp(-(Frequenzen - om).^2 ./ (4*MD.Gamma));
    %Gaborfenster     = 2 * exp(-Gamma^2 * (Frequenzen - om).^2);
    Transformed(i,X) = ifft( Y .* Gaborfenster);
end;
