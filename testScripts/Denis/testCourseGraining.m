clear all

%Signal
fs=250; %Hz
Ts = 1/fs; %sec
Nt = 50000; %number of points
t = [1:Nt].'*Ts; %time in sec

x = randn(Nt,1); %gaussian noise
SignalName = 'WhiteNoise';
% x = zeros(Nt,1);
% for iF = 1:10;
%     x = x + sin(2*pi*((iF*10+1)*t));
% end
% SignalName = 'MultiFreauencyHarmonicOscillators';

%Scaling
scales = [1:10:51];
Nsc = length(scales);
y{1} = x;
time{1} = t;
for iSc = 2:Nsc;
    sc2 = floor(scales(iSc)/2);
    time{iSc} = [];
    y{iSc} = [];
    for iP = (sc2+1):scales(iSc):Nt-sc2;
        time{iSc} = [time{iSc};t(iP)];
        %Coarse graining
        y{iSc} = [y{iSc}; mean( x(iP-sc2:iP+sc2) ) ];
         %Downsampling
%         y{iSc} = [y{iSc}; x(iP) ];
         %Resampling
%         y(iP,iSc) = resample(x(iP),floor(Nt/scales(iSc)),Nt);
        
    end
end



%Power spectra
f = 1:100;
for iSc = 1:Nsc;
    %For white noise:
    P(:,iSc)= pwelch(y{iSc},hamming(length(y{iSc})/10),[],f,fs/scales(iSc));
%     %For harmonic oscillators:
%   P(:,iSc) = pwelch(y{iSc},length(y{iSc}),[],[1:100],fs/scales(iSc));
end


%Plotting
%cmap = colormap(jet(Nsc));
figure('Name',SignalName);
for iSc = 1:Nsc;
    subplot(Nsc,2,2*(iSc-1)+1)
    plot(time{iSc},y{iSc}); 
    axis tight;hold off;
    ylabel(['Scale: ',num2str(scales(iSc)*Ts),' sec'])
    subplot(Nsc,2,2*(iSc-1)+2)
    plot(f,P(:,iSc));
    axis([f(1) f(end) 0 max(P(:,iSc))]);hold off;
end
subplot(Nsc,2,1)
hold on;
title('Time series')
hold off;
subplot(Nsc,2,2)
hold on;
title('Power spectra')
hold off;
subplot(Nsc,2,2*Nsc-1)
hold on;
xlabel('Time (sec)');
hold off;
subplot(Nsc,2,2*Nsc)
hold on;
xlabel('Frequency (Hz)');
hold off;
