%clear all

%Create 10 trials of gaussian noise:
x= [abs(randn(10001,1)), rand(10001,1)];%randn(1001,10);

%Set some parameters:
%Obligatory:
cfg.fs=1000; %sampling frequency in Herz
cfg.method = 'trial'; %other options are 'trialTime', and 'ensemble'
cfg.measures= {'DFA','SP','VGR','MSEall','AC'};%(alternatively 'all')

%Optional:
cfg.time = [0:0.001:10].';%[0:0.001:1].'; %time vector
cfg.normal = 'zscore'; %recommended if method='trial'

scales = 1:10;
%Structures for measures' parameters
%MSEall:
cfg.MSEall.m = 2;
cfg.MSEall.r = 0.5;
cfg.MSEall.scales = scales;
%DFA:
cfg.DFA.order = 2;
cfg.DFA.scales = scales(scales>cfg.DFA.order+1);
%VGR: 
cfg.VGR.scales = scales;
%AC: 
cfg.AC.scales = scales;
%cfg.AC.norm = 1;

%Run configuration function:
[cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x);
%Run measures' calculation function:
[C, cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x);
