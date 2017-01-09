clear all

%Create 10 trials of gaussian noise:
x= [abs(randn(10001,1)), rand(10001,1)];%randn(1001,10);

%Set some parameters:
%Obligatory:
cfg.fs=1000; %sampling frequency in Herz
cfg.method = 'trial'; %other options are 'trialTime', and 'ensemble'
cfg.measures= {'PDF','CDF'};%{'NE','SE','MSE','DFA'}; %(alternatively 'all')

%Optional:
cfg.time = [0:0.001:10].';%[0:0.001:1].'; %time vector
cfg.normal = 'none'; %recommended if method='trial'

%Structures for measures' parameters
%NE:
%cfg.NE.Nbins = sqrt(1000);
%cfg.NE.approach = 'biased';
%cfg.NE.base = exp(1);
%cfg.NE.norm = 1;
%SE:
%cfg.SE.m = 2;
%cfg.SE.r = 0.5;
%MSE:
%cfg.MSE.m = 2;
%cfg.MSE.r = 0.5;
%cfg.MSE.scales = [1:20];
%DFA:
%cfg.DFA.order = 2;
%cfg.DFA.scales = [1:20];

cfg.PDF.method = 'kernelisohist';
cfg.CDF.method = 'interpol';
%Run configuration function:
[cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x);
%Run measures' calculation function:
[C, cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x);
% %--------------------------------------------------------------------------
% %In case that:
% cfg.method = 'ensemble';%'trialTime';%'ensemble';
% cfg.timeCalc = cfg.time([1:50:951]); %points of calculation in time, e.g., centers of time windows of calculation
% cfg.winLen = 0.1; %length of time window in secs (time)
% % cfg.timeCalc = cfg.time(501);
% % cfg.winLen = 1;
% 
% %Run configuration function:
% [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, x] = DPtimeUnivarPrepare(cfg,x);
% 
% %Run measures' calculation function:
% [C, cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,x);

%--------------------------------------------------------------------------
% %For statistics:
% cfg.stats.method = 'timeshuffling';%'phaseshuffling','phaserandom'
% cfg.stats.alpha = 0.01;
% cfg.stats.tail = 2;
% cfg.stats.Nperm = 100;
% cfg.stats.pointStatMethod = '';%'z','t'
% cfg.stats.corrMultComp = 'FDR'; %'BONF';
% 
% %Run configuration function:
% [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x)
% 
% %Run measures' calculation function:
% [C, cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures,thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x);