function [C, cfg, statsRes] = DPtimeUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas, x)

%%This function prepares data and configuration struture for DPstateCoupling
%

%Inputs: 
%--------------------------------------------------------------------------
%-cfg: configuration structure with fields
%   -fs: sampling frequency in herz, positive real scalar
%   -time: time vector in secs, real valued vector, default=[0:N-1]/fs
%--------------------------------------------------------------------------
%   -method: one of 
%            1. 'trial', for estimation per trial
%            2. 'trialTime', for estimation per trial with a sliding time 
%               window
%            3. 'ensemble', for estimation across trials, either pointwise
%            or with a time window
%--------------------------------------------------------------------------
%   Optional, but recommended if method is 'trial':
%   -normal: 'zscore', 'meanCntr', 'linear', or 'none' string, default='none' 
%           if normal == 'linear', there should also be a vector of 2 real
%           numbers, normVal for the minimum and maximum value, 
%           default, normVal = [-1 1]
%--------------------------------------------------------------------------
%Optionally for 'trialTime' or 'ensemble'
%   -timeCalc: a subset of timeP (in secs), with the time points where the calculation
%       will be performed, default=timeP
%   -winLen: a positive number for the time length of the time window of
%   calculation, default = 1/fs 
% leading to a window length of Nwin points
%
%when timeCalc is a donwsampled version of time or timeP, then
%   -upsample: 'yes' or 'no', in order to interpolate back to
%       the initial time vector
%
%either after a calculation with timeCalc=timeP or after upsampling and
%           if smoothing is desired
%   -smoothWinfun: a function handle to a MATLAB window function 
%            in case of smoothing after calculation, default: @hanning
%   -smoothWinlen: a positive number for the time length of the time window of
%               smoothing, default = min(Tc)/4, i.e. a quarter of the
%               central period or main time scale of the slower signal
%--------------------------------------------------------------------------
%   -measures: the coupling measures that should be calculated 
%              among the available ones:
%              1.'NE' for naive entropy calculation (based on 1D histogram)
%              2.'SE' for sample entropy
%              3.'MSE' for multiscale entropy 
%              4.'MSEall' for all variants of MSE, except for CMSE
%              5.'CMSE' for composite sample entropy
%              6.'DFA' for detrended fluctuation analysis
%              7.'SP' for spectral measures: DOF (Degress Of
%                 Freedom), PS (power spectral slope in a log-log graph), 
%                 logP (logarithm of spectral power).
%              8.'VGR' for variogram and its log-log slope
%              9.'PDF' for probability density function
%              10.'CDF' for cumulative probability density function
%              11.'AC' for autocorrelation function
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
%       -downsample: takes values in {0,1}, 1 for using matlab function "resample" and
%                    0 for the traditional coarse-graining, default=0;
%       -normScale: takes values in {0,1}, if normScale=1, data are 
%                   normalized for each scale separately, default=0;
%       -slopeInds: the indices of the scale points to be used for the
%                   calculation of the slope,  vector of positive integers 
%                   monotonously increasing, default: 1:length(scales)
%                   empty matrix to skip calculation of slope and return
%                   zero
%   -MSEall: optional structure of inputs related to MSE
%       -m: number points of the pattern, positive integer scalar
%           default=2
%       -r: similarity threshold, real, positive, scalar in the interval (0,1)
%           default=0.5
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/50)
%       -slopeInds: the indices of the scale points to be used for the
%                   calculation of the slope,  vector of positive integers 
%                   monotonously increasing, default: 1:length(scales)
%                   empty matrix to skip calculation of slope and return
%                   zero
%   -CMSE: optional structure of inputs related to MSE
%       -m: number points of the pattern, positive integer scalar
%           default=2
%       -r: similarity threshold, real, positive, scalar in the interval (0,1)
%           default=0.5
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/50)
%       -normScale: takes values in {0,1}, if normScale=1, data are 
%                   normalized for each scale separately, default=0;
%       -slopeInds: the indices of the scale points to be used for the
%                   calculation of the slope,  vector of positive integers 
%                   monotonously increasing, default: 1:length(scales)
%                   empty matrix to skip calculation of slope and return
%                   zero
%   -DFA: optional structure of inputs related to DFA
%       -order: order of the DFA, positive, scalar, integer, default=2
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: order+2:floor(N/10)
%       -slopeInds: the indices of the scale points to be used for the
%                   calculation of the slope, vector of positive integers
%                   monotonously increasing, default: 1:length(scales),
%                   empty matrix to skip calculation of slope and return
%                   zero
%   -SP:  optional structure of inputs related to SP
%       -NFFT: the number of frequency bins for the FFT, postive integer
%              scalar and a power of 2, in the interval (2, nextpow2(N)],
%              default=nextpow2(N)
%       -winfun: a function handle to a window function, default: (@N)boxcar(N)
%       -slopeInds: the indices of the frequency points to be used for the
%                   calculation of the slope, vector of positive integers 
%                   monotonously increasing, default: 1:length(scales)
%                   empty matrix to skip calculation of slope and return
%                   zero
%   -VGR: optional structure of inputs related to VGR
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/10)
%       -slopeInds: the indices of the scale points to be used for the
%                   calculation of the slope,  vector of positive integers 
%                   monotonously increasing, default: 1:length(scales)
%                   empty matrix to skip calculation of slope and return
%                   zero
%   -PDF: optional structure of inputs related to PDF
%       -method: one of the following strings:
%                -'isodist': equidistant histogram (or other histogram
%                where bin centers in signal's space are specified)
%                -'isohist': iso-probability histogram (based on quantiles)
%                -'kernelisodist': like 'isodist' but with kernel smoothing
%                -'isohist': like 'isohist' but with kernel smoothing
%                default is 'isodist'
%       -bins: either the number of bins (a positive integer scalar >=3) 
%              or a vector of bin centers for methods '(kernel)isodist'
%                 (a vector of monotonously increasing real values)
%              or a vector of quantiles for methods '(kernel)isohist'
%                 (a vector of monotonously increasing real values 
%                  in the interval [0 1])
%              default bins = sqrt(N)
%      -norm: equal to 1 in order to normalize pdf to have a total area of
%             1, otherwise equal to 0, default is 1
%      -kernel: in case that a kernel smoothing method is selected
%               kernel can be a string or a function handle to a kernel
%               function (see help of ksdensity.m)
%               default = 'normal'
%      -width:  in case that a kernel smoothing method is selected 
%               it is the bandwidth of the kernel (a positive real number
%               <= max(x)-min(x))
%               default = smallest bin width in signal's space
%   -CDF: optional structure of inputs related to CDF
%       -method: one of the following strings:
%                -'interpol': equidistant histogram (or other histogram
%                where bin centers in signal's space are specified), linear
%                interpolation of the empirical cumulative distribution
%                function is performed (see help of ecdf.m)
%                -'quantl': iso-probability histogram (based on quantiles)
%                -'kernelinterpol': like 'interpol' but with kernel smoothing
%                -'kernelquantl': like 'quantl' but with additional kernel smoothing
%                default is 'quantl'
%       -bins: either the number of bins (a positive integer scalar >=3) 
%              or a vector of bin edges for methods '(kernel)interpol'
%                 (a vector of monotonously increasing real values)
%              or a vector of quantiles for methods '(kernel)quantl'
%                 (a vector of monotonously increasing real values 
%                  in the interval [0 1])
%              default bins = sqrt(N)
%      -kernel: in case that a kernel smoothing method is selected
%               kernel can be a string or a function handle to a kernel
%               function (see help of ksdensity)
%               default = 'normal'
%      -width:  in case that a kernel smoothing method is selected 
%               it is the bandwidth of the kernel (a positive real number
%               <= max(x)-min(x))
%               default = smallest bin width in signal's space
%   -AC: optional structure of inputs related to AC
%       -scales: scales' vector in number of points, vector of positive,
%                integers, default: 1:floor(N/10)
%       -slopeInds: the indices of the scale points to be used for the
%                   calculation of the slope,  vector of positive integers 
%                   monotonously increasing, default: 1:length(scales)
%                   empty matrix to skip calculation of slope and return
%                   zero
%      -norm: equal to 1 in order to normalize AC to have a maximum of 1
%             for scale 0, otherwise equal to 0, default is 1
%--------------------------------------------------------------------------
%Optionally if statistics are to be made:
%   -stats: a structure with fields:
%       -method:a string with the surrogate test method that is to be applied:
%                 -"timeshuffling": random time shuffling 
%                 -"phaseshuffling": random phase shuffling in Fourier
%                 space
%                 -"phaserandom": randomized phase in Fourier space
%                default: "timeshuffling"
%       -surrfun: a handle to a function of the form
%                 xSurr = surrfun(x), where
%                   -x: the data series columnwise matrix, 
%                that creates a surrogate data set given the method above
%       -alpha: alpha value, positive real number<=0.5, default: 0.05
%       -tail: 1 or 2 for 1- or 2-talied tests, default: 1
%--------------------------------------------------------------------------
%      Optionally, in the case of multiple comparisons correction
%      and if pointStat exists:
%       -corrMultComp: 'BONF' or 'FDR', default:'' for no correction
%--------------------------------------------------------------------------
%       -Nperm: number of permutations, positive integer,
%               equal or greater than the default value:
%               tail*ceil(1/alpha)
%       -pointStatMethod: a string defining the point statistic to be used,
%                         among 't' (for t-values),
%                               'z' for z-score normalization, or '' (none),
%                         default = '' 
%      Optionally,
%       -multiStatfun: a structure with fields:
%                      -function: a handle to a function of the form
%                                 multiStat = multiStatfun(C,Csurr,pointStat,params) where
%                                 -C: the result of the original data
%                                 -Csurr: the results of the surrogate data
%                                 -pointStat: a structure of the form of C
%                                             with values of the point statistic
%                                 -multiStat: a structure of the form of C
%                                             with values a multivariate statistic
%                     -params: parameters for that function
%--------------------------------------------------------------------------": random phase shuffling in Fourier
%                 space
%                 -"phaserandom": randomized phase in Fourier space
%                default: "timeshuffling"
%       -alpha: alpha value, positive real number<=0.5, default: 0.01
%       -tail: 1 or 2 for 1- or 2-talied tests, default: 1
%       -Nperm: number of permutations, positive integer,
%               equal or greater than the default value: tail*ceil(1/alpha)+1
%       -pointStatMethod: a string defining the point statistic to be used,
%                         among 't' (for t-values),
%                               'z' for z-score normalization, or '' (none),
%                         default = '' 
%      Optionally,
%       -multiStatfun: a structure with fields:
%                      -function: a handle to a function of the form
%                                 multiStat = multiStatfun(C,Csurr,pointStat,params) where
%                                 -C: the result of the original data
%                                 -Csurr: the results of the surrogate data
%                                 -pointStat: a structure of the form of C
%                                             with values of the point statistic
%                                 -multiStat: a structure of the form of C
%                                             with values a multivariate statistic
%                     -params: parameters for that function
%      Optionally, in the case of multiple comparisons correction
%      and if pointStat exists:
%       -corrMultComp: 'BONF' or 'FDR', default:'' for no correction
%--------------------------------------------------------------------------


% %  -MeasNames: the cell with the names of all (sub)measures in the results'
% %             structure as constructed here:
% 
% % -NmeasPmeas: the vector of the numbers of submeasures per measure 
% %              in the results' structure as constructed here:
% -NoutsPmeas: the number of outputs per calculation for each measure
% measures' categories
%  -method: 1,2,or 3 depending on whether cfg.method =  'trial', 
%           'trialTime', or 'ensemble' respectively 
%  -measInds: the indexes of the measures to be calculated
%  -Nmeasures: how many measures are to be calculated
%  -thisfunCommands: cell of the command strings for the calculation of
%                    the selected measures
%  -x: input data processed (e.g., z-scored), same size as x in the input 



%  Outputs:

%  -C: results' structure with fields with the names of the calculated measures
%       of size depending on the method:
%       1. vector 1 x M for 'method='trial'
%       2. matrix cfg.Ncalc x M for 'method='trialTime'
%       3. vector Nout x 1  for method='ensemble'
%     if method == 'trial' or 'trialTime', then mean values across trials are
%     given in a substructure C.trialMean as:
%       1. scalar for 'method='trial'
%       2. vector cfg.Nout x 1 for 'trialTime'
%  -cfg: configuration structure corrected and complemented
%  -statsRes: results cell structures for the surrogate statistics' test,
%             same structure as C, with fieds
%             -poinStat: values of the statistic per data point
%             -multiStat: values of multivariate statistic
%             -p: p values
%             -pUncorr: p values without correction for multiple
%                       comparisons
%             -sign: flags of value 1 or 0 for significant or non
%                    significant results

% tic

%Main calculation
switch method
    
    
    case 1 %'trial'
        
        %The dimensions of the output of each measure
        cfg.Nout(1) = cfg.Ntr;
        
        %The time vector of the output
        cfg.timeOut = cfg.time;
        
        %Initialize measure vectors:
        %For each measure we will calculate...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                Nouts = NoutsPmeas{measInds(iM)}(iSubM);
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ntr,Nouts);
                if cfg.Ntr>1
                    C.trialMean.(MeasNames{measInds(iM),iSubM})=zeros(1,Nouts);
                end
            end
        end
        
        
        %The calculation loop:
        
        %...for each trial...
        for iTr=1:cfg.Ntr;
            
            %disp(['...calculating trial',num2str(iT),'/',num2str(cfg.Ntr),'...'])
            
            %...get the time series of this trial
            thisX = x(:,iTr);
            
            %...calculate the selected measures for this trial
            for iC = 1:Ncomnds;
                try
                    thisfunCommands{iC}
                    eval(thisfunCommands{iC});
                catch
                    keyboard;
                end
            end
            
        end
        
        
        %...calculate the mean of all trials:
        if cfg.Ntr>1
            %For each measure we calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...calculate the mean value (excluding possible nan
                    %values)
                    %C.trialMean.(MeasNames{measInds(iM),iSubM}) = nanmean( C.(MeasNames{measInds(iM),iSubM}) );
                    C.trialMean.(MeasNames{measInds(iM),iSubM}) = mean( C.(MeasNames{measInds(iM),iSubM}) );
                end
            end
        end
        
        
        
    case 2 %'trialTime'
        
        %The dimensions of the output of each measure
        cfg.Nout(1) = cfg.Ncalc;
        cfg.Nout(3) = cfg.Ntr;
        
        %The time vector of the output
        cfg.timeOut = cfg.timeCalc;
        
        %Initialize measure vectors:
        %For each measure we will calculate...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                Nouts = NoutsPmeas{measInds(iM)}(iSubM);
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,Nouts,cfg.Ntr);
            end
        end
        
        
        
        %The calculation loop:
        
        %Initialize time window centers and edges for this u
        cfg.timeointWinEdgs = zeros(cfg.Ncalc,2);
        cfg.timeWinEdgs = zeros(cfg.Ncalc,2);
        
        
        %For each calculation/window step...
        for iT=1:cfg.Ncalc;
            
            %...calculate indexes for:
            cfg.timeointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
            cfg.timeointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.N);%ending window point
            cfg.timeWinEdgs(iT,:) = cfg.time(cfg.timeointWinEdgs(iT,:)); %same points in time...
            
            
            %...for each trial...
            for iTr=1:cfg.Ntr;
                
                %...get the time series of this trial and time window
                thisX = x( cfg.timeointWinEdgs(iT,1):cfg.timeointWinEdgs(iT,2),iTr );
                N=length(thisX);
                
                %...calculate the selected measures for this trial and time point
                for iC = 1:Ncomnds;
                    try
                        thisfunCommands{iC}
                        eval(thisfunCommands{iC});
                    catch
                        keyboard
                    end
                end
                
                
            end
            
        end
        
        %For each measure we calculated...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure ...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                
                if cfg.Ntr>1
                    %...calculate the mean value (excluding possible nan
                    %values)...
                    %C.trialMean.(MeasNames{measInds(iM),iSubM}) = squeeze(nanmean( C.(MeasNames{measInds(iM),iSubM}), 3 ));
                    C.trialMean.(MeasNames{measInds(iM),iSubM}) = squeeze(mean( C.(MeasNames{measInds(iM),iSubM}), 3 ));
                end
                %...squeeze the trials' matrix
                C.(MeasNames{measInds(iM),iSubM})=squeeze(C.(MeasNames{measInds(iM),iSubM}));
                
            end
        end
        %Interpolate through upsampling if required:
        if strcmpi(cfg.upsample,'yes')
            oldC = C.trialMean; %temporary copy
            %For each measure we calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure ...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...upsample from Ncalc points to N:
                    C.trialMean.(MeasNames{measInds(iM),iSubM}) = resample( oldC.(MeasNames{measInds(iM),iSubM}), cfg.N, cfg.Ncalc );
                end
            end
            cfg.Nout(1) = cfg.N;
            cfg.timeOut = cfg.time;
        end
        clear oldC;
        
        %Smooth through convolution with a window function if required:
        if isfield(cfg,'smoothWin')
            %For each measure we calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure ...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...smooth:
                    C.trialMean.(MeasNames{measInds(iM),iSubM}) = conv2( C.trialMean.(MeasNames{measInds(iM),iSubM}), cfg.smoothWin, 'same' );
                end
            end
        end
        
        
    otherwise %'ensemble'
        
        %The dimensions of the output of each measure
        cfg.Nout(1) = cfg.Ncalc;
        
        %The time vector of the output
        cfg.timeOut = cfg.timeCalc;
        
        %Initialize measure vectors:
        %For each measure we will calculate...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                Nouts = NoutsPmeas{measInds(iM)}(iSubM);
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,Nouts);
            end
        end
        
        %Initialize time window centers and edges
        cfg.timeointWinEdgs = zeros(cfg.Ncalc,2);
        cfg.timeWinEdgs = zeros(cfg.Ncalc,2);
        
        
        %For each calculation/window step...
        for iT=1:cfg.Ncalc;
            
            %...calculate indexes for:
            cfg.timeointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
            cfg.timeointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.N);%ending window point
            cfg.timeWinEdgs(iT,:) = cfg.time(cfg.timeointWinEdgs(iT,:)); %same points in time...
            
            %...get the time series of this trial and time window
            thisX = x( cfg.timeointWinEdgs(iT,1):cfg.timeointWinEdgs(iT,2) , :) ;
            N=size(thisX,1);
            
            %...calculate the selected measures for this trial and time point
            for iC = 1:Ncomnds;
                try
                    thisfunCommands{iC}
                    eval(thisfunCommands{iC});
                catch
                    keyboard
                end
            end
            
        end
        
        
        
        %Interpolate through upsampling if required:
        if strcmpi(cfg.upsample,'yes')
            oldC = C; %temporary copy
            %For each measure we calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure ...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...upsample from Ncalc points to N:
                    C.(MeasNames{measInds(iM),iSubM}) = resample( oldC.(MeasNames{measInds(iM),iSubM}), cfg.N, cfg.Ncalc ); 
                end
            end
            cfg.Nout(1) = cfg.N;
            cfg.timeOut = cfg.time;
        end
        
        %Smooth through convolution with a window function if required:
        if isfield(cfg,'smoothWin')
            %For each measure we calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure ...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...smooth:
                    C.(MeasNames{measInds(iM),iSubM}) = conv2( C.(MeasNames{measInds(iM),iSubM}), cfg.smoothWin, 'same' );
                end
            end
        end
        
%         toc
        
end


%For the statistics..


%In case surrogate statistics are to be made...
if ~isempty(cfg.stats)
    
    
    %.unload parameters...
    pointStatMethod=cfg.stats.pointStatMethod;
    multiStatfun=cfg.stats.multiStatfun;
    alpha=cfg.stats.alpha;
    tail=cfg.stats.tail;
    corrMultComp=cfg.stats.corrMultComp;
    
    %...initialize a cell to store the measures structures for the
    %surrogates...
    Csurr=cell(1,cfg.stats.Nperm);
    
    %...and a configurations' structure, identical to cfg, except for
    %the surrogate statistics...
    cfgSurr = cfg;
    cfgSurr.stats=[];
    
    %...for of the surrogate time series...
    for iST=1:cfg.stats.Nperm;
        
        %...create it...
        xSurr = cfg.stats.surrfun(x,1,cfg.Ntr); 

        %...and calculate the measures with exactly the same setting
        %except for the surrogate statistics...
        disp(['...calculating measures for surrogate data set ',num2str(iST),'/',num2str(cfg.stats.Nperm),'...'])
        [Csurr{iST}, dummy, dummy2] = DPtimeUnivar(cfgSurr, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, NoutsPmeas,xSurr); 
        
    end
    
    
    %Initialize statistics results structure
    statsRes = struct();
    
    %...initialize measure vectors:
    %For each measure we have calculated...
    for iM = 1:Nmeasures;
        %...and for each of the sumbmeasures of this measure...
        for iSubM = 1:NmeasPmeas(measInds(iM));
                            
            disp(['...calculating statistics for measure ',MeasNames{measInds(iM),iSubM},'...'])
            
            Nouts = NoutsPmeas{measInds(iM)}(iSubM);
            
            %...select method...
            switch method
                
                case 1 %trial
                    
                    
                    %...create a temporary matrix of the measure of the surrogate dataset...
                    thisCsurr =  zeros(cfg.Ntr,Nouts,cfg.stats.Nperm);
                    for iST=1:cfg.stats.Nperm;
                        %...for each permutation...
                        thisCsurr(:,:,iST) = Csurr{iST}.(MeasNames{measInds(iM),iSubM});
                    end
                    thisCsurr = squeeze(thisCsurr);
                    
                    %...calculate the specific statistic for each trial...
                    [pointStat,  dummy] = DPcalcStat(C.(MeasNames{measInds(iM),iSubM}),thisCsurr,pointStatMethod, multiStatfun);
                    clear thisCsurr;
                    
                    %...calculate significance for each trial...
                    [statsRes.(MeasNames{measInds(iM),iSubM}).p,...
                     statsRes.(MeasNames{measInds(iM),iSubM}).pUnCorr,...
                     statsRes.(MeasNames{measInds(iM),iSubM}).sign ] = ...
                     ...
                     DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                    
                    %...store statistics except for the ones of the surrogates...
                    %pointStat = rmfield(pointStat,'surrVal');
                    statsRes.(MeasNames{measInds(iM),iSubM}).stat = pointStat;
                    
                    
                    %...if trials are more than one do the same for the trials' mean...
                    if cfg.Ntr>1
                        
                        %...create a temporary matrix of the measure of the surrogate dataset...
                        thisCmeanSurr =  zeros(Nouts,cfg.stats.Nperm);
                        for iST=1:cfg.stats.Nperm;
                            thisCmeanSurr(:,iST) = Csurr{iST}.trialMean.(MeasNames{measInds(iM),iSubM});
                        end
                        thisCmeanSurr = squeeze(thisCmeanSurr);
                        
                        %...calculate the specific statistic...
                        [pointStat,  dummy] = DPcalcStat(C.trialMean.(MeasNames{measInds(iM),iSubM}),thisCmeanSurr,pointStatMethod, multiStatfun);
                        clear thisCmeanSurr;
                        
                        %...and calculate significance...
                        [ statsRes.(MeasNames{measInds(iM),iSubM}).pMean,...
                          statsRes.(MeasNames{measInds(iM),iSubM}).pMeanUnCorr,...
                          statsRes.(MeasNames{measInds(iM),iSubM}).signMean ] = ...
                          ...
                          DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                        
                        %...store statistics except for the ones of the surrogates...
                        %pointStat = rmfield(pointStat,'surrVal');
                        statsRes.(MeasNames{measInds(iM),iSubM}).meanStat = pointStat;
                        
                    end
                    
                    
                    
                case 2 %trialTime
                    
                    %...create a temporary matrix of the measure of the surrogate dataset...
                    thisCsurr =  zeros(cfg.Ncalc,Nouts,cfg.Ntr,cfg.stats.Nperm);
                    for iST=1:cfg.stats.Nperm;
                        %...for each trial...
                        thisCsurr(:,:,:,iST) = Csurr{iST}.(MeasNames{measInds(iM),iSubM});
                    end
                    thisCsurr = squeeze(thisCsurr);
                    
                    %...calculate the specific statistic for each trial...
                    [pointStat,  multiStat] = DPcalcStat(C.(MeasNames{measInds(iM),iSubM}),thisCsurr,pointStatMethod, multiStatfun);
                    clear thisCsurr;
                    
                    %...calculate significance for each trial for each time point...
                    [statsRes.(MeasNames{measInds(iM),iSubM}).p,...
                     statsRes.(MeasNames{measInds(iM),iSubM}).pUnCorr,...
                     statsRes.(MeasNames{measInds(iM),iSubM}).sign ] = ...
                     ...
                     DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                    
                    %...store statistics except for the ones of the surrogates...
                    %pointStat = rmfield(pointStat,'surrVal');
                    statsRes.(MeasNames{measInds(iM),iSubM}).pointStat = pointStat;
                    
                    
                    if ~isempty(multiStat)
                        
                        %...calculate significance for each trial for multivariate statistic (all times)...
                        [statsRes.(MeasNames{measInds(iM),iSubM}).pMulti,...
                         statsRes.(MeasNames{measInds(iM),iSubM}).pMultiUnCorr,...
                         statsRes.(MeasNames{measInds(iM),iSubM}).signMulti ] = ...
                         ...
                         DPcalcSignif(multiStat,alpha,tail,corrMultComp);
                        
                        %...store statistics except for the ones of the surrogates...
                        %multiStat = rmfield(multiStat,'surrVal');
                        statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).multiStat = multiStat;
                    end
                    
                    %...if trials are more than one do the same for the trials' mean...
                    if cfg.Ntr>1
                        
                        %...create a temporary matrix of the measure of the surrogate dataset...
                        thisCmeanSurr =  zeros(cfg.Ncalc,Nouts,cfg.stats.Nperm);
                        for iST=1:cfg.stats.Nperm;
                            thisCmeanSurr(:,:,iST) = Csurr{iST}.trialMean.(MeasNames{measInds(iM),iSubM});
                        end
                        thisCmeanSurr = squeeze(thisCmeanSurr);
                        
                        %...calculate the specific statistic...
                        [pointStat,  multiStat] = DPcalcStat(C.trialMean.(MeasNames{measInds(iM),iSubM}),thisCmeanSurr,pointStatMethod, multiStatfun);
                        clear thisCmeanSurr;
                        
                        %...and calculate significance for each time point....
                        [ statsRes.(MeasNames{measInds(iM),iSubM}).pMean,...
                          statsRes.(MeasNames{measInds(iM),iSubM}).pMeanUnCorr,...
                          statsRes.(MeasNames{measInds(iM),iSubM}).signMean ] = ...
                          ...
                          DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                        
                        %...store statistics except for the ones of the surrogates...
                        %pointStat = rmfield(pointStat,'surrVal');
                        statsRes.(MeasNames{measInds(iM),iSubM}).meanPointStat = pointStat;
                        
                        if ~isempty(multiStat)
                            
                            %...calculate significance for multivariate statistic (all times)...
                            [ statsRes.(MeasNames{measInds(iM),iSubM}).pMultiMean,...
                              statsRes.(MeasNames{measInds(iM),iSubM}).pMultiMeanUnCorr,...
                              statsRes.(MeasNames{measInds(iM),iSubM}).signMultiMean ] = ...
                              ...
                              DPcalcSignif(multiStat,alpha,tail,corrMultComp);
                            
                            %...store statistics except for the ones of the surrogates...
                            %multiStat = rmfield(multiStat,'surrVal');
                            statsRes.(MeasNames{measInds(iM),iSubM}).meanMultiStat = multiStat;
                            
                        end
                    end
                    
                    
                    
                otherwise
                    
                    %...create a temporary matrix of the measure of the surrogate dataset...
                    thisCsurr =  zeros(cfg.Ncalc,Nouts,cfg.stats.Nperm);
                    for iST=1:cfg.stats.Nperm;
                        thisCsurr(:,:,iST) = Csurr{iST}.(MeasNames{measInds(iM),iSubM});
                    end
                    thisCsurr = squeeze(thisCsurr);
                    
                    %...calculate the specific statistic...
                    [pointStat,  multiStat] = ...
                        DPcalcStat(C.(MeasNames{measInds(iM),iSubM}),thisCsurr,pointStatMethod, multiStatfun);
                    clear thisCsurr;
                    
                    %...calculate significance for each time point...
                    [statsRes.(MeasNames{measInds(iM),iSubM}).p,...
                     statsRes.(MeasNames{measInds(iM),iSubM}).pUnCorr,...
                     statsRes.(MeasNames{measInds(iM),iSubM}).sign ] = ...
                     ...
                     DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                    
                    %...store statistics except for the ones of the surrogates...
                    %pointStat = rmfield(pointStat,'surrVal');
                    statsRes.(MeasNames{measInds(iM),iSubM}).pointStat = pointStat;
                    
                    if ~isempty(multiStat)
                        %...calculate significance for multivariate statistic (all times)...
                        [statsRes.(MeasNames{measInds(iM),iSubM}).pMulti,...
                         statsRes.(MeasNames{measInds(iM),iSubM}).pMultiUnCorr,...
                         statsRes.(MeasNames{measInds(iM),iSubM}).signMulti ] = ...
                         ...
                         DPcalcSignif(multiStat,alpha,tail,corrMultComp);
                        
                        %...store statistics except for the ones of the surrogates...
                        %multiStat = rmfield(multiStat,'surrVal');
                        statsRes.(MeasNames{measInds(iM),iSubM}).multiStat = multiStat;
                    end
                    
                    
            end %method selection
 
        end %submeasure
    end %measure
    
else
    statsRes=[];
end

%disp('DONE!')



