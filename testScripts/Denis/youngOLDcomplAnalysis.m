function youngOLDcomplAnalysis(filename)

%This script will calculate the selected univariate complexity measures in  
%time domain, with the selected method for ONE ONLY dataset comprising a
%specific number of trials and EEG channels.
%One could loop this script across subjects and conditions...

%Toolbox locations: /Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox
%Put those into MATLAB path:
% DataHandling
%io: input output between VA and MATLAB
%Statistics (if statistics are to be made)
%Univariate Toolbox
%Utilities

% clear all;
% clc;


%A. Load VA exported Data into MATLAB
%(!export from VA as binary with format IEEE_FLOAT_32 and multiplexed!)
%Inputs:
%filename = 'mdo13117_1RestEC10s250Hz.bin';% full path, filename and extension of the VA data file
%Optional:
%channels: a vector of channel numbers [1 2 4 etc], or string 'all'
%trials: a vector of trials/segmetns numbers [1 2 4 etc], or string 'all'
[data hdr mrk] = DPreadBV(filename);

%data has and should have the format time x channels x trials/segments
% 
[N Nch Ntr] = size(data);
%or
% N=hdr.SegmentDataPoints;
% Ntr = hdr.nTrials;
% Nch = hdr.NumberOfChannels;


%B. Parameters of configutation structure:

%1. Obligatory:
cfg.fs=hdr.Fs; %Sampling frequency in Hz
cfg.method = 'trial';%(1 value per trial, no time window) 
                         %other: 'trialTime' (a sliding window within trial), 
                         %or 'ensemble' (sliding window across trials)

%2. Optional but recommended:

%   Optional, but recommended if method is 'trial':
%   -normal: 'zscore', 'meanCenter', 'linear', or 'none' string, default='zscore' 
%           if normal == 'linear', there should also be a vector of 2 real
%           numbers, normVal for the minimum and maximum value, 
%           default, normVal = [-1 1]
% cfg.normal='zscore';

cfg.time = [0:N-1]/hdr.Fs;% in seconds, default: cfg.time = [0 N-1]/cfg.fs;

%3. In case that 'trialTime', or 'ensemble' is selected:
%The time points (sliding window centers) where the calculation will be
%performed. It has to be a subset of cfg.time. default: cfg.timeCalc =
%cfg.time, for pointwise calculation.
% cfg.timeCalc = cfg.time(ceil(N/2)); 
% cfg.winLen = 10-1/cfg.fs; %the window length in seconds, it will always have an odd number of points

%4. More optional parameters:

%In case of a sliding window:
%when timeCalc is a donwsampled version of time or timeCut, then
%   -upsample: 'yes' or 'no', in order to interpolate back to
%       the initial time vector
%
%either after a calculation with timeCalc=timeCut or after upsampling and
%           if smoothing is desired
%   -smoothWinfun: a function handle to a MATLAB window function 
%            in case of smoothing after calculation, default: @hanning
%   -smoothWinlen: a positive number for the time length of the time window of
%               smoothing, default = min(Tc)/4, i.e. a quarter of the
%               period of the slower signal

%5. Measures and their parameters:
%--------------------------------------------------------------------------
%   -measures: the coupling measures that should be calculated 
%              among the available ones:
%              1.'NE' for naive entropy calculation (based on 1D histogram)
%              2.'SE' for sample entropy
%              3.'MSE' for multiscale entropy 
%              4.'CMSE' for composite sample entropy
%              5.'DFA' for detrended fluctuation analysis
%              6.'SP' for spectral measures: DOF (Degress Of
%                 Freedom), PS (power spectral slope in a log-log graph), 
%                 logP (logarithm of spectral power).
%              7.'VGR' for variogram and its log-log slope
%
%   -NE: optional structure of inputs related to NE
%       -Nbins: the number of bins for the histogram, positive integer scalar,
%               default = round(sqrt(N))
%       -approach: The method used, one of the following ones:
%                  'biased': The biased estimate (default)
%                  'mmse': The minimum mean square error estimate
%                  'unbiased' : The unbiased estimate
%       -base: the base of the logarithm, positive, real, scalar, default e
%       -norm: 1, for normalization of the output in the interval [0 1] (default)
%              0, otherwise, in which case value and units depend on base
%   -SE: optional structure of inputs related to SE
%       -m: number points of the pattern, positive integer scalar
%           default=2
%       -r: similarity threshold, real, positive, scalar in the interval (0,1)
%           default=0.5
%   -MSE: optional structure of inputs related to MSE
%       -m: number points of the pattern, positive integer scalar
%           default=2
%       -r: similarity threshold, real, positive, scalar in the interval (0,1)
%           default=0.5
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/50)
%       -normScale: takes values in {0,1}, if normScale=1, data are 
%                   normalized for each scale separately, default=0;
%   -CMSE: optional structure of inputs related to MSE
%       -m: number points of the pattern, positive integer scalar
%           default=2
%       -r: similarity threshold, real, positive, scalar in the interval (0,1)
%           default=0.5
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/50)
%       -normScale: takes values in {0,1}, if normScale=1, data are 
%                   normalized for each scale separately, default=0;
%   -DFA: optional structure of inputs related to DFA
%       -order: order of the DFA, positive, scalar, integer, default=2
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: order+2:floor(N/10)
%   -SP:  optional structure of inputs related to SP
%       -NFFT: the number of frequency bins for the FFT, postive integer
%              scalar and a power of 2, in the interval (2, nextpow2(N)],
%              default=nextpow2(N)
%       -winfun: a function handle to a window function, default: (@N)boxcar(N)
%   -VGR: optional structure of inputs related to VGR
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/10)

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
%C. Optionally if statistics are required:
% cfg.stats.method = 'timeshuffling'; %Method of surrogates: 
%                                 Oher options: % -"timeshuffling": random time shuffling 
%                                                 -"phaseshuffling": random phase shuffling in Fourier space
%                                                 -"phaserandom": randomized phase in Fourier space
%                                                 default: "timeshuffling"
%cfg.stats.tail =1; %default =1, otherwise, =2 for two tailed test
%Other options:
%      Optionally, in the case of multiple comparisons correction
%      and if pointStat exists:
%       -corrMultComp: 'BONF' or 'FDR', default:'' for no correction
%       -Nperm: number of permutations, positive integer,
%               equal or greater than the default value:
%               tail*ceil(1/alpha)
%       -pointStatMethod: a string defining the point statistic to be used,
%                         among 't' (for t-values),
%                               'z' for z-score normalization, or '' (none),
%                         default = '' 

%
%      Optionally,
%       -multiStatfun: a structure with fields:
%                      -function: a handle to a function of the form
%                                 multiStat = multiStatfun(PC,PCsurr,pointStat,params) where
%                                 -PC: the result of the original data
%                                 -PCsurr: the results of the surrogate data
%                                 -pointStat: a structure of the form of PC
%                                             with values of the point statistic
%                                 -multiStat: a structure of the form of PC
%                                             with values a multivariate statistic
%                     -params: parameters for that function


%D. calculation loop for all electrodes


%Choose a subset of channels if you wish...
Channels = [1:Nch];% %60;
Nch=length(Channels);

%Make one calculation to get the size of the results' matrices for each
%measure
%[cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, x] =  DPtimeUnivarPrepare(cfg,x);
%Inputs:
x = squeeze(data(:,1,:));
cfgTemp  = cfg;
cfgTemp.stats = [];
%[cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x)
[cfgTemp, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfgTemp,x);
%cfg.stats = [];

%Initialize measure matrices for method 'trial'
for iM = 1:Nmeasures;
    %...and for each of the sumbmeasures of this measure...
    for iSubM = 1:NmeasPmeas(measInds(iM));
        Nouts = NoutsPmeas{measInds(iM)}(iSubM);
        %The result is given in a structure in the form:
        %C.MeasureName  = time (window) (x scales -for MSE and DFA) x trial
        %matrix
        %RES.(MeasNames{measInds(iM),iSubM})=zeros(cfgTemp.Ncalc,Nouts,Ntr,Nch);
        RES.(MeasNames{measInds(iM),iSubM})=zeros(Ntr,Nouts,Nch);
        %RES.(MeasNames{measInds(iM),iSubM})=zeros(cfgTemp.Ncalc,Nouts,Nch);
        if cfgTemp.Ntr>1
            %and for the trial's mean (unless 'ensemble' is chosen):
            %C.trialMean.MeasureName  = time (window) (x scales -for MSE and DFA- x)
            %matrix
            REStrialMean.(MeasNames{measInds(iM),iSubM})=zeros(Nouts,Nch);
            %REStrialMean.(MeasNames{measInds(iM),iSubM})=zeros(cfgTemp.Ncalc,Nouts,Nch);
        end
    end
end
% statsRes = cell(1,Nch);
clear cfgTemp;

%for all channels
for ii=1:Nch;
    
   % tic
    
    iC=Channels(ii);
    
    iC
    
    %[cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] =  DPtimeUnivarPrepare(cfg,x);
    %Inputs:
    x = squeeze(data(:,iC,:));
    [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x);
    
    %[C, cfg, statsRes, MeasNames, NmeasPmeas, NoutsPmeas] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,x)
    [C, cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x);
    
    
    %Unpack results and put them in a common structure for all channels
    for iM = 1:Nmeasures;
        %...and for each of the sumbmeasures of this measure...
        for iSubM = 1:NmeasPmeas(measInds(iM));

            %For method 'trialTime', the result is given in a shape of the form of a matrix:
            %C.MeasureName  = time (window) x 1 (or x scales -for (C)MSE, DFA and VGR- or frequencies for SP) x trial:
            %RES.(MeasNames{measInds(iM),iSubM})(:,:,:,ii)=C.(MeasNames{measInds(iM),iSubM});
            
            %Instead, for method 'trial' we have 1 (or scales -for (C)MSE, DFA
            %and VGR- or frequencies for SP) x trial matrix:
            RES.(MeasNames{measInds(iM),iSubM})(:,:,ii)=C.(MeasNames{measInds(iM),iSubM});
            %RES.(MeasNames{measInds(iM),iSubM})(:,:,ii)=C.(MeasNames{measInds(iM),iSubM});
            %...and for method 'ensemble' we have:
            %Instead, for method 'trial' we have time (window) x1 (or x scales -for (C)MSE, DFA and VGR- or frequencies for SP):
            
            %RES.(MeasNames{measInds(iM),iSubM})(:,:,ii)=C.(MeasNames{measInds(iM),iSubM});
            
            if cfg.Ntr>1
                %and for the trial's mean for 'trialTime':
                %REStrialMean.(MeasNames{measInds(iM),iSubM})(:,:,ii)=C.trialMean.(MeasNames{measInds(iM),iSubM});
                
                %and for 'trial':
                REStrialMean.(MeasNames{measInds(iM),iSubM})(:,ii)=C.trialMean.(MeasNames{measInds(iM),iSubM});
                

            end
        end
    end
    
    %Accordingly for the results of the statistics, one finds inside:
    %p values,
    %p_UnCorr for NO correction for multiple comparisons
    %sign: a matric of 1s (or 0s) for the values that are (not) significant
    %in the same shape as the result
    % and some of the statistics of the surrogate data set
%toc
end

%Squeeze measure matrices
for iM = 1:Nmeasures;
    %...and for each of the sumbmeasures of this measure...
    for iSubM = 1:NmeasPmeas(measInds(iM));

        %The result is given in a structure in the form:
        %C.MeasureName  = trial x scales -for MSE and DFA) 
        %matrix
        RES.(MeasNames{measInds(iM),iSubM})=squeeze(permute(RES.(MeasNames{measInds(iM),iSubM}),[2,3,1]));
        %RES.(MeasNames{measInds(iM),iSubM})=squeeze(RES.(MeasNames{measInds(iM),iSubM}));
        if cfg.Ntr>1
            %and for the trial's mean (unless 'ensemble' is chosen):
            %C.trialMean.MeasureName  = time (window) (x scales -for MSE and DFA- x)
            %matrix
            REStrialMean.(MeasNames{measInds(iM),iSubM})=squeeze(REStrialMean.(MeasNames{measInds(iM),iSubM}));
        end
    end
    
end
datasetName = hdr.DataFile(1:end-4);
save(['C:\data\DATA_Aging_Complexity\Results\RES',cfg.method,'_',datasetName,'.mat'],'cfg','RES','REStrialMean','method', 'measInds', 'Nmeasures', 'thisfunCommands', 'Ncomnds', 'MeasNames', 'NmeasPmeas', 'NoutsPmeas')%,'statsRes'



%E. Optionally one can save these results now as Vision Analyzer files in
%order to load the into VA and visualize them, e.g. topographies...


