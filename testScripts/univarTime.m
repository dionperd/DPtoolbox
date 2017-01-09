%Toolbox locations: /Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox
%Put those into MATLAB path:
%io: input output between VA and MATLAB
%Statistics (if statistics are to be made)
%Univariate Toolbox
%Utilities
%Data Handling

clear all;
clc;


%A. Load VA exported Data into MATLAB
%(!export from VA as binary with format IEEE_FLOAT_32 and multiplexed!)
%Inputs:
filename = './DPdata.eeg';% full path, filename and extension of the VA data file
%Optional:
%channels: a vector of channel numbers [1 2 4 etc], or string 'all'
%trials: a vector of trials/segmetns numbers [1 2 4 etc], or string 'all'
[data hdr mrk mrkPtype] = DPreadBV(filename);
%data has and should have the format time x channels x trials/segments



%B. Parameters of configutation structure:

%1. Obligatory:
cfg.fs=hdr.Fs; %Sampling frequency in Hz
cfg.method = 'trialTime';%this means a sliding window within trial,
                         %other: 'trial' (1 value per trial, no time window), 'ensemble' (sliding window across trials)

%2. Optional but recommended:
%cfg.measures = {'MSE', 'NE'}; %default: 'all' 
cfg.time = [1:563]*0.001;% in seconds, default: cfg.time = [0 N-1]/cfg.fs;

%3. In case that 'trialTime', or 'ensemble' is selected:
%The time points (sliding window centers) where the calculation will be
%performed. It has to be a subset of cfg.time. default: cfg.timeCalc =
%cfg.time, for pointwise calculation.
cfg.timeCalc = cfg.time([100, 250, 400]); 
cfg.winLen = 0.5; %the window length in seconds, it will always have an odd number of points

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
%              1.'NE' for naive shannon entropy calculation (based on 1D histogram)
%              2.'SE' for sample entropy
%              3.'MSE' for multiscale sample entropy 
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
%       -NFFT: the number of frequency bins for the FFT, positive integer
%              scalar and a power of 2, in the interval (2, nextpow2(N)],
%              default=nextpow2(N)
%   -VGR: optional structure of inputs related to VGR
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/10)

%C. Optionally if statistics are required:
cfg.stats.method = 'timeshuffling'; %Method of surrogates: 
%                                 Oher options: % -"timeshuffling": random time shuffling 
%                                                 -"phaseshuffling": random phase shuffling in Fourier space
%                                                 -"phaserandom": randomized phase in Fourier space
%                                                 default: "timeshuffling"
%cfg.stats.alpha = 0.05; %default = 0.05
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



%D. calculation loop for all frequencies -not CFC- and electrode pairs



Channels = [1 30 60]; %60;
Nch=length(Channels);

%for all channels
for ii=1:Nch;
    
    iC=Channels(ii);
    
    %[cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] =  DPtimeUnivarPrepare(cfg,x);
    %Inputs:
    x = squeeze(data(:,iC,:));
    [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x);
    
    %[C, cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x)
    [C, cfg, statsRes{ii}] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x);
    
    
    %Unpack results and put them in a common structure for all channels
    for iM = 1:Nmeasures;
        %...and for each of the sumbmeasures of this measure...
        for iSubM = 1:NmeasPmeas(measInds(iM));
            Nouts = NoutsPmeas{measInds(iM)}(iSubM);
            %For method 'trialTime', the result is given in a shape of the form of a matrix:
            %C.MeasureName  = time (window) x 1 (or x scales -for (C)MSE, DFA and VGR- or frequencies for SP) x trial:
            RES.(MeasNames{measInds(iM),iSubM})(:,:,:,iC)=C.(MeasNames{measInds(iM),iSubM});
            
            %Instead, for method 'trial' we have 1 (or scales -for (C)MSE, DFA
            %and VGR- or frequencies for SP) x trial matrix:
     
            %...and for method 'ensemble' we have:
            %Instead, for method 'trial' we have time (window) x1 (or x scales -for (C)MSE, DFA and VGR- or frequencies for SP):
            
            %RES.(MeasNames{measInds(iM),iSubM})(:,:,iC)=C.(MeasNames{measInds(iM),iSubM});
            
            if cfg.Ntr>1
                %and for the trial's mean for 'trialTime':
                REStrialMean.(MeasNames{measInds(iM),iSubM})(:,:,iC)=C.trialMean.(MeasNames{measInds(iM),iSubM});
                
                %and for 'trial':
                %REStrialMean.(MeasNames{measInds(iM),iSubM})(:,iC)=C.trialMean.(MeasNames{measInds(iM),iSubM});
                
            end
        end
    end
    
    %Accordingly for the results of the statistics, one finds inside:
    %p values,
    %p_UnCorr for NO correction for multiple comparisons
    %sign: a matric of 1s (or 0s) for the values that are (not) significant
    %in the same shape as the result
    % and some of the statistics of the surrogate data set

end

%Squeeze measure matrices
for iM = 1:Nmeasures;
    %...and for each of the sumbmeasures of this measure...
    for iSubM = 1:NmeasPmeas(measInds(iM));
        Nouts = NoutsPmeas{measInds(iM)}(iSubM);
        %The result is given in a structure in the form:
        %C.MeasureName  = time (window) (x scales -for MSE and DFA) x trial
        %matrix
        RES.(MeasNames{measInds(iM),iSubM})=squeeze(RES.(MeasNames{measInds(iM),iSubM}));
        if cfg.Ntr>1
            %and for the trial's mean (unless 'ensemble' is chosen):
            %C.trialMean.MeasureName  = time (window) (x scales -for MSE and DFA- x)
            %matrix
            REStrialMean.(MeasNames{measInds(iM),iSubM})=squeeze(REStrialMean.(MeasNames{measInds(iM),iSubM}));
        end
    end
end

save(['univarTime',cfg.method,'.mat'],'RES','REStrialMean','statsRes')






