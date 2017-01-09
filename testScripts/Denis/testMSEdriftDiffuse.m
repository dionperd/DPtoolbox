function testMSEdriftDiffuse

t1 = tic;



%Parameters of complexity measures calculation
cfg.fs=250; %Sampling frequency in Hz
dt=1/cfg.fs;
cfg.method = 'trial';
cfg.time = [dt:dt:12]/hdr.Fs;% in seconds, default: cfg.time = [0 N-1]/cfg.fs;
Nt = length(cfg.time);
cfg.measures = {'MSE','DFA','SP','VGR'}; %'SE','CMSE'
% cfg.SE.m=2;
% cfg.SE.r=0.5;
cfg.MSE.m=2;
cfg.MSE.r=0.5;
cfg.MSE.scales = [1:50];
cfg.MSE.normScale = 0; 
%Better to calculate CMSE than MSE... more robust...
% cfg.CMSE.m=2;
% cfg.CMSE.r=0.5;
% cfg.CMSE.scales = [1:50];
% cfg.CMSE.normScale = 0; %normalize std per scale or not???
%Better to calculate DFA than power spectra slopes
%The variogram is also a nice option and much faster than DFA...
cfg.DFA.order=2;
cfg.DFA.scales=[4:50];
cfg.VGR.scales = [1:50];
cfg.SP.winfun=@(N)hamming(N);
cfg.SP.NFFT=2^nextpow2(N);

%Parameters of integrations:
cfg.method = 'Euler';
cfg.systemType = 'nonauto';
cfg.s_noise=0;
plotting = -10;
dt=0.001;
t = [dt:dt:12].';
Nt = length(t);
fs = 1000; %Hz
Noise = randn(Nt,1);

Conds = {'Linear 1D', 'Non linear 2D oscillator'};
Groups = {'s=100', 's=10', 's=1'};
T = [0.01 0.1 1];
sigma = 1./T;
Ns = length(sigma);
for iS = 1:Ns;

    %Linear 1D stochastic system
    pL.t=t;
    pL.T=T(iS);
    pL.sigma=sigma(iS);
    pL.Noise = Noise;
    System=@(t,x,p)Linear1D(t,x,p);
    x0=0;
    disp('Simulating linear 1D system...')
    tic
    [ t, xL, pout, dxdt, dxdt_ns, cfg] = integration_NoMex(System,pL,t,x0,cfg,plotting+2);
    toc
    if any(isnan(xL))
        error('Unstable integration!')
    end
    %Downsample to fs Hz
    xL = resample(xL,fs*t(end),Nt);
    
    %Oscillatory 2D stochastic system
    pO.t=t;
    pO.T=T(iS);
    pO.mu = 1;
    pO.sigma=sigma(iS);
    pO.Noise = [Noise, randn(size(Noise))];
    System=@(t,x,p)Oscillatory2D(t,x,p);
    x0=[0; 0];
    disp('Simulating nonlinear 2D oscillatory system...')
    tic
    [ t, xO, pout, dxdt, dxdt_ns, cfg] = integration_NoMex(System,pO,t,x0,cfg,plotting+3);
    toc
    if any(isnan(xO))
        error('Unstable integration!')
    end
    xO = xO(:,1);
    %Downsample to fs Hz
    xO = resample(xO,fs*t(end),Nt);
    
    x = [xL, xO];
    t = t(1):(t(end)-t(1))/(Nt-1):t(end);
    
    %Get rid of transients and edge effects;
    ind = (t>=1) & (t<=(t(end)-1));
    x=x(ind,:);
    t=t(ind);
    Nt=length(t);
    
    %Normalize
    x=zscore(x);
    %x=x-repmat(mean(x),[N,1]);
    
    %Calculate measures
    [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x);
    [C(iS), cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x);

    
    
    %Plotting
    
    
    for iS=1:2;
        figure(ht);
        subplot(5,1,iS)
        title(label{iS})
        plot(t,x(:,iS),'Color',Colors{iS},'linewidth',1);
        hold off;
        
        figure(hf);
        subplot(1,2,1)
        hold on;
        plot(log(f),log(P(:,iS)),'Color',Colors{iS},'linewidth',1);
        
        subplot(1,2,2)
        hold on;
        plot(f,P(:,iS),'Color',Colors{iS});
        
    end
    figure(ht)
    subplot(5,1,5);
    hold on;
    xlabel('Time (sec)');
    hold off;
    
    figure(hf)
    subplot(1,2,1)
    hold on;
    legend(label);
    xlabel('log(f) (log(Hz))');
    grid on;
    %axis equal
    hold off;
    subplot(1,2,2)
    hold on;
    xlabel('f (Hz)');
    grid on;
    hold off;
    
end
toc(t1)



function [dxdt, pout] = Linear1D(t,x,p)
pout = p.sigma*p.Noise(t==p.t);
dxdt =  -x/p.T + pout;

function [dxdt, pout] = Oscillatory2D(t,x,p)
pout = p.sigma*p.Noise(t==p.t,:);
dxdt(1,1) =  ( x(1) + x(2)-x(1)^3 )/p.mu/p.T  + pout(1);
dxdt(2,1) =  -p.mu*( x(1) - 0.5*x(2) )/p.T + pout(2) ;

