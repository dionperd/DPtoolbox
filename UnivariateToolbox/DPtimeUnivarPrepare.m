function [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, NoutsPmeas, x] = DPtimeUnivarPrepare(cfg,x)

%This function prepares data and configuration struture for DPtimeUnivar
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



%  Outputs:
%  -cfg: configuration structure corrected and complemented
%  -method: 1,2,or 3 depending on whether cfg.method =  'trial', 
%           'trialTime', or 'ensemble' respectively 
%  -measInds: the indexes of the measures to be calculated
%  -Nmeasures: how many measures are to be calculated
%  -thisfunCommands: cell of the command strings for the calculation of
%                    the selected measures
% %  -MeasNames: the cell with the names of all (sub)measures in the results'
% %             structure as constructed here:

% 
% % -NmeasPmeas: the vector of the numbers of submeasures per measure 
% %              in the results' structure as constructed here:
% -NoutsPmeas: the number of outputs per calculation for each measure
%  -x: input data processed (e.g., z-scored), same size as x in the input 




%Data validation, some unpacking and calculation of parameters and commands
funcName='DPtimeUnivarPrepare';


xSIZE = size(x);
xSIZElen = length(xSIZE);
varName='x';
testX = {@(x)isnumeric(x),...
         @(x)isreal(x),...
         @(x,xSIZElen)xSIZElen<=3 };

param={{},{},{xSIZElen}};
mode=['e','e','e'];
execfun={{},{},{}};
default=nan;
[RESULT, x] = DPvalidateData(x,testX,param,mode,execfun,default,varName,funcName);
x=squeeze(x);
N = xSIZE(1); %Number of time points
cfg.N =N;
Ntr=xSIZE(2);
cfg.Ntr =Ntr;



varName='cfg';
testCFG = {@(cfg)isstruct(cfg) };
param = {{}};
mode=['e'];
execfun={{}};
default=nan;
[RESULT, cfg] = DPvalidateData(cfg,testCFG,param,mode,execfun,default,varName,funcName);

if isfield(cfg,'fs')
    varName='fs';
    testFS = {@(fs)isnumeric(fs),...
              @(fs)isreal(fs),...
              @(fs)isscalar(fs),...
              @(fs)( fs>0 ) };
    param={{},{},{},{}};
    mode=['e','e','e','e'];
    execfun={{},{},{},{}};
    default=nan;
    [RESULT, cfg.fs] = DPvalidateData(cfg.fs,testFS,param,mode,execfun,default,varName,funcName);
    cfg.Ts=1/cfg.fs; %sampling time in secs
else
    error('Structure field cfg.fs is missing.')
end

%   -fc: central frequency of the signals in herz, 2 (or 3) element vector of
%        positive real values
% if isfield(cfg,'fc')
%     varName='fc';
%     testFC = {@(fc)isnumeric(fc),...
%               @(fc)isreal(fc),...
%               @(fc)all( fc>0 ),...
%               @(fc)isvector(fc),...
%               @(fc,MultiVar)( (MultiVar==0) && (length(fc)==2) ) || ( (MultiVar==1) && (length(fc)==3) ) };
%     param={{},{},{},{},{cfg.MultiVar}};
%     mode=['e','e','e','e','e'];
%     execfun={{},{},{},{},{}};
%     default=nan;
%     [RESULT, cfg.fc] = DPvalidateData(cfg.fc,testFC,param,mode,execfun,default,varName,funcName);
%     cfg.Tc=1./cfg.fc; %central period in secs
%     cfg.TcMax = max(cfg.Tc); %and for the slower signal in case of CFC
%     cfg.NtcMax = round(cfg.TcMax/cfg.Ts); %and in time points (samples)
% else
%     error('Structure field cfg.fc is missing.')
% end


if isfield(cfg,'time')
    varName='time';
    testTIME = { @(time)isnumeric(time),...
                 @(time)isreal(time),...
                 @(time)isvector(time),...
                 @(time)all(diff(time)>0),...
                 @(time,N)length(time)<=N};
    param={{},{},{},{},{N}};
    mode=['e','e','w','e','e'];
    execfun={{},{},@(time)time(:),{},{}};
    default=([0:N-1]*cfg.Ts).';
    [RESULT, cfg.time] = DPvalidateData(cfg.time,testTIME,param,mode,execfun,default,varName,funcName);
else
    %cprintf('Magenta','WARNING: Structure field cfg.time is missing. Setting default [0:N-1]/fs.\n')
    fprintf('WARNING: Structure field cfg.time is missing. Setting default [0:N-1]/fs.\n')
    cfg.time = ([0:N-1]*cfg.Ts).';
end
cfg.timeLen = cfg.time(end)-cfg.time(1);

if isfield(cfg,'method')
    varName='method';
    testMETHOD = { @(method)ischar(method),...
                   @(method)isvector(method),...
                   @(method)any(strcmpi(method,{'trial','trialTime','ensemble'})) };
    param={{},{},{}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default=nan;
    [RESULT, cfg.method] = DPvalidateData(cfg.method,testMETHOD,param,mode,execfun,default,varName,funcName);
else
    error('Structure field cfg.method is missing.')
end
method = 1*strcmpi(cfg.method,'trial') + 2*strcmpi(cfg.method,'trialTime') + 3*strcmpi(cfg.method,'ensemble');


%Normalization only if method is 'trial' and there is some specific method
%requested
if strcmpi(cfg.method,'trial')
    
    if isfield(cfg,'normal')
        varName='normal';
        testNORMAL = { @(normal)ischar(normal),...
                       @(normal)isvector(normal),...
                       @(normal)any(strcmpi(normal,{'none','zscore','meanCntr','linear'})) };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='none';
        [RESULT, cfg.normal] = DPvalidateData(cfg.normal,testNORMAL,param,mode,execfun,default,varName,funcName);
        
        if strcmpi(cfg.normal,'linear')
            if isfield(cfg,'normVal')
                varName='normVal';
                testNORMVAL = { @(normVal)isnumeric(normVal),...
                    @(normVal)isvector(normVal),...
                    @(normVal)length(normVal)==2,...
                    @(normVal)normVal(2)>normVal(1) };
                param={{},{},{},{}};
                mode=['e','e','e','e'];
                execfun={{},{},{},{}};
                default=[-1 1];
                [RESULT, cfg.normVal] = DPvalidateData(cfg.normVal,testNORMVAL,param,mode,execfun,default,varName,funcName);
            else
                cfg.normVal = [-1 1];
            end
        end
    else
        cfg.normal = 'zscore';
    end
    
    if strcmpi(cfg.normal,'zscore')
        
        [x,cfg.xM,cfg.xSTD] = zscore(x);
        
    elseif strcmpi(cfg.normal,'meanCntr')
        
        cfg.xM = mean(x);
        x = x-repmat(cfg.xM,N,1);
        
    elseif strcmpi(cfg.normal,'linear')
        
        cfg.minX = min(x);
        cfg.maxX = max(x);
        range = cfg.maxX - cfg.minX;
        x = normVal(1) + repmat((diff(normalVal)./range),N,1) ./ ( x - repmat(cfg.minX,N,1) );
        
    end
    
    %DP: too complicated and not well suited for the measures we calculate for
    %the moment here
    % if strcmpi(cfg.normal,'zscore')
    %     if strcmpi(cfg.method, 'ensemble')
    %        [x,cfg.xM,cfg.xSTD] = zscore(x(:));
    %     else
    %         [x,cfg.xM,cfg.xSTD] = zscore(x);
    %     end
    % elseif strcmpi(cfg.normal,'meanCntr')
    %     if strcmpi(cfg.method, 'ensemble')
    %         cfg.xM = mean(x);
    %         x = x-cfg.xM ;
    %     else
    %         cfg.xM = mean(x);
    %         x = x-repmat(cfg.xM,N,1);
    %     end
    % elseif strcmpi(cfg.normal,'linear')
    %     if strcmpi(cfg.method, 'ensemble')
    %         cfg.minX = min(x(:));
    %         cfg.maxX = max(x(:));
    %         range = cfg.maxX - cfg.minX;
    %         x = normVal(1) + repmat((diff(normalVal)./range),N,1) ./ ( x - cfg.minX );
    %     else
    %         cfg.minX = min(x);
    %         cfg.maxX = max(x);
    %         range = cfg.maxX - cfg.minX;
    %         x = normVal(1) + repmat((diff(normalVal)./range),N,1) ./ ( x - repmat(cfg.minX,N,1) );
    %     end
    % end
    
else
    cfg.normal='none';
end


if (method~=1)  
    
    if isfield(cfg,'timeCalc')
        varName='timeCalc';
        testTIMECALC = { @(timeCalc)isnumeric(timeCalc),...
                         @(timeCalc)isreal(timeCalc),...
                         @(timeCalc)isvector(timeCalc),...
                         @(timeCalc)all(diff(timeCalc)>0),...
                         @(timeCalc,N)length(timeCalc)<=N,...
                         @(timeCalc,time) all( ismember(timeCalc,time) )};
        param={{},{},{},{},{cfg.N},{cfg.time}};
        mode=['e','e','w','e','e','e'];
        execfun={{},{},@(timeCalc)timeCalc(:),{},{},{}};
        default=cfg.time;
        [RESULT, cfg.timeCalc] = DPvalidateData(cfg.timeCalc,testTIMECALC,param,mode,execfun,default,varName,funcName);
    else
        %cprintf('Magenta','WARNING: Structure field cfg.timeCalc is missing. Setting default cfg.timeCalc = cfg.timeCut.\n')
        fprintf('WARNING: Structure field cfg.timeCalc is missing. Setting default cfg.timeCalc = cfg.timeCut.\n')
        cfg.timeCalc = cfg.time;
    end
    cfg.Ncalc = length(cfg.timeCalc);
    
    %Initialize the time indices of calculation to be equal with all times:
    cfg.timeCalcIND = 1:cfg.N;
    if ~isequal(cfg.timeCalc,cfg.time)
        if isfield(cfg,'upsample')
            varName='upsample';
            testUPSAMPLE = { @(upsample)ischar(upsample),...
                             @(upsample)isvector(upsample),...
                             @(upsample)any(strcmpi(upsample,{'yes','no'})) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='no';
            [RESULT, cfg.upsample] = DPvalidateData(cfg.upsample,testUPSAMPLE,param,mode,execfun,default,varName,funcName);
        else
            cfg.upsample='no';
            %cprintf('Magenta','WARNING: Structure field cfg.upsample is missing. Setting default ''no''.\n')
            fprintf('WARNING: Structure field cfg.upsample is missing. Setting default ''no''.\n')
        end
        %Update cfg.timeCalcIND
        cfg.timeCalcIND = cfg.timeCalcIND( ismember(cfg.time,cfg.timeCalc) );
    else
        cfg.upsample='no';
    end
    
    if isfield(cfg,'winLen')
        varName='winLen';
        testWINLEN = { @(winLen)isnumeric(winLen),...
                       @(winLen)isreal(winLen),...
                       @(winLen)isscalar(winLen),...
                       @(winLen,timeLen)(winLen>=0)&&(winLen<=timeLen)};
        param={{},{},{},{cfg.timeLen}};
        mode=['e','e','e','e'];
        execfun={{},{},{},{}};
        default=cfg.Ts;
        [RESULT, cfg.winLen] = DPvalidateData(cfg.winLen,testWINLEN,param,mode,execfun,default,varName,funcName);
    else
        %cprintf('Magenta','WARNING: Structure field cfg.winLen is missing. Setting cfg.winLen=Ts=%f.\n',cfg.Ts)
        fprintf('WARNING: Structure field cfg.winLen is missing. Setting cfg.winLen=Ts=%f.\n',cfg.Ts)
        cfg.winLen = cfg.Ts;
    end
    %Find window length in points
    cfg.Nwin = ceil(cfg.winLen / cfg.Ts);
    %We need the window to have odd number of points
    if mod(cfg.Nwin,2)~=1
        cfg.Nwin=cfg.Nwin+1;
    end
    cfg.Nwin2 = floor(cfg.Nwin/2); %a useful constant, the half window 
    %Update window length in time
    cfg.winLen = cfg.Nwin*cfg.Ts;
    
    if isfield(cfg,'smoothWinfun')&&( isequal(cfg.timeCalc,cfg.time) || strcmpi(cfg.upsample,'yes') )
        varName='smoothWinfun';
        testSMOOTHWINFUN = { @(smoothWinfun)isa(smoothWinfun,'function_handle') };
        param={{}};
        mode=['e'];
        execfun={{}};
        default=@hanning;
        [RESULT, cfg.smoothWinfun] = DPvalidateData(cfg.smoothWinfun,testSMOOTHWINFUN,param,mode,execfun,default,varName,funcName);

        if isfield(cfg,'smoothWinlen')
            varName='smoothWinlen';
            testSMOOTHWINLEN = { @(smoothWinlen)isnumeric(smoothWinlen),...
                                 @(smoothWinlen)isreal(smoothWinlen),...
                                 @(smoothWinlen)isscalar(smoothWinlen),...
                                 @(smoothWinlen,Ts,timeLen)(smoothWinlen>Ts)&&(smoothWinlen<timeLen)};
            param={{},{},{},{cfg.Ts,cfg.timeLen}};
            mode=['e','e','e','e'];
            execfun={{},{},{},{}};
            default=min(cfg.Tc)/4;
            [RESULT, cfg.smoothWinlen] = DPvalidateData(cfg.smoothWinlen,testSMOOTHWINLEN,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.smoothWinlen is missing. Setting cfg.smoothWinlen=min(cfg.Tc)/4=%f.\n',min(cfg.Tc)/4)
            fprintf('WARNING: Structure field cfg.smoothWinlen is missing. Setting cfg.smoothWinlen=min(cfg.Tc)/4=%f.\n',min(cfg.Tc)/4)
            cfg.smoothWinlen = min(cfg.Tc)/4;
        end
        %Find smoothing window length in points
        cfg.NsmoothWin = ceil(cfg.smoothWinlen / cfg.Ts);
        %We need the window to have odd number of points
        if mod(cfg.NsmoothWin,2)~=1
            cfg.NsmoothWin=cfg.NsmoothWin+1;
        end
        cfg.NsmoothWin2 = floor(cfg.NsmoothWin/2); %a useful constant, the half window
        %Update window length in time
        cfg.smoothWinlen = cfg.NsmoothWin*cfg.Ts;
        
        %Calculate the smoothing window
        cfg.smoothWin=cfg.smoothWinfun(cfg.Nwin);
        %Normalize... 
        cfg.smoothWin=cfg.smoothWin/sum(cfg.smoothWin);
        %... and vectorize it as well
        cfg.smoothWin = cfg.smoothWin(:);
    end
end


if isfield(cfg,'measures')
    varName='measures';
    testMEASURES = { @(measures)iscellstr(measures)||ischar(measures),...
                     @(measures)~isempty(intersect(measures,{'all','NE','SE','MSE','MSEall','CMSE','DFA','SP','PDF','CDF','AC'})) };
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    default='all';
    [RESULT, cfg.measures] = DPvalidateData(cfg.measures,testMEASURES,param,mode,execfun,default,varName,funcName);
    
else
    %cprintf('Magenta','WARNING: Structure field cfg.measures is missing. Replacing with ''all''.\n')
    fprintf('WARNING: Structure field cfg.measures is missing. Replacing with ''all''.\n')
    cfg.measures='all';
end


%  -MeasNames: the cell with the names of all (sub)measures in the results'
%             structure as constructed here:
MeasNames={'NE','','','','','','','','','';... %NE, naive (Shannon) entropy
           'SE','','','','','','','','','';... %SE, sample entropy
           'MSE','STD','G','','','','','','','';... %MSE, multiscale sample entropy and std
           'MSEc','RMSEc','STDc','MSEr','RMSEr','STDr','Gc','RGc','Gr','RGr';... %MSE, multiscale sample entropy and std, all variants
           'CMSE','CSTD','CG','','','','','','','';... %CMSE, composite multiscale sample entropy and std
           'H','logF','','','','','','','','';... %DFA, detrended fluctuation analysis
           'P','DOF','logP','PS','','','','','','';...% Spectral analysis
           'logV','VS','','','','','','','',''; %Variogram
            'PDF','PDFxi'  ,'','','','','','','','';
            'CDF','CDFxi'  ,'','','','','','','','';
            'R','AG'  ,'','','','','','','',''; %Autocorrelation and its log slope
            };

% -NmeasPmeas: the vector of the numbers of submeasures per measure 
%              in the results' structure as constructed here:
%             NE SE MSE MSEall CMSE DFA SP VGR PDF CDF AC
NmeasPmeas = [ 1  1   3   10     3   2   4   2   2  2  2 ]; %Number of 'measures per measures' categories


%Initialization:
% -NoutsPmeas: the cell of the numbers of output points per submeasure in
%              the results' structure, as constructed here:
%             NE SE    MSE     MSEall      CMSE     DFA       SP      V       PDF      CDF    AC
NoutsPmeas = {1; 1;  [1 1 1]; ones(1,10); [1 1 1]; [1 1]; ones(1,4); [1 1] ;  [1 1]  ; [1 1]; [1 1]};
NoutsMAX=1;

%Unpack the measures that will be calculated
if strcmpi(cfg.measures,'all')
    measures = ones(1,11);
else
    measures = [any(strcmpi(cfg.measures,'NE')),...
                any(strcmpi(cfg.measures,'SE')),...
                any(strcmpi(cfg.measures,'MSE')),...
                any(strcmpi(cfg.measures,'MSEall')),...
                any(strcmpi(cfg.measures,'CMSE')),...                
                any(strcmpi(cfg.measures,'DFA')),...
                any(strcmpi(cfg.measures,'SP')),...
                any(strcmpi(cfg.measures,'VGR')),... 
                any(strcmpi(cfg.measures,'PDF')),... 
                any(strcmpi(cfg.measures,'CDF')),...
                any(strcmpi(cfg.measures,'AC')),...
                ];
end


%Calculate command strings
if (method==1) %trial
    
    funCommands  = {'[C.(''NE'')(iTr), estimation, nbias,sigma] = DPcalcNE(thisX,cfg.NE.Nbins,cfg.NE.approach,cfg.NE.base,cfg.NE.normInp,cfg.NE.norm);', ...
                    'C.(''SE'')(iTr) = DPcalcSampEn(thisX,cfg.SE.m,cfg.SE.r,cfg.N,cfg.SE.norm);', ...
                    '[C.(''MSE'')(iTr,:), C.(''STD'')(iTr,:),C.(''G'')(iTr,:)] = DPcalcMSE(thisX,cfg.MSE.m,cfg.MSE.r,cfg.MSE.scales,cfg.N,1,cfg.MSE.Nsc,,cfg.MSE.downsample, cfg.MSE.normScale,cfg.MSE.slopeInds);', ...
                    '[C.(''MSEc'')(iTr,:), C.(''RMSEc'')(iTr,:), C.(''STDc'')(iTr,:), C.(''MSEr'')(iTr,:), C.(''RMSEr'')(iTr,:), C.(''STDr'')(iTr,:),C.(''Gc'')(iTr,:),C.(''RGc'')(iTr,:),C.(''Gr'')(iTr,:),C.(''RGr'')(iTr,:)] = DPcalcMSEall(thisX,cfg.MSEall.m,cfg.MSEall.r,cfg.MSEall.scales,cfg.N,1,cfg.MSEall.Nsc,cfg.MSEall.slopeInds);', ...
                    '[C.(''CMSE'')(iTr,:), C.(''CSTD'')(iTr,:),C.(''CG'')(iTr,:)]] = DPcalcCMSE(thisX,cfg.CMSE.m,cfg.CMSE.r,cfg.CMSE.scales,cfg.N,1,cfg.CMSE.Nsc,cfg.cfg.CMSE.normScale,cfg.CMSE.slopeInds);', ...                    
                    '[C.(''H'')(iTr,:), C.(''logF'')(iTr,:) ] = DPcalcDFA(thisX,cfg.DFA.scales,cfg.DFA.logScales,cfg.DFA.order,cfg.N,cfg.DFA.Nsc,cfg.DFA.slopeInds);',...
                    '[C.(''P'')(iTr,:), C.(''DOF'')(iTr), C.(''logP'')(iTr,:), C.(''PS'')(iTr,:) ] = DPcalcSP(thisX,cfg.SP.logf,cfg.SP.winfun,cfg.SP.NFFT,cfg.SP.Nf,cfg.N,cfg.SP.slopeInds);',...
                    '[C.(''logV'')(iTr,:), C.(''VS'')(iTr,:) ] = DPcalcVarGram(thisX,cfg.VGR.scales,cfg.VGR.logScales,cfg.VGR.Nsc,cfg.VGR.slopeInds);',...
                    '[C.(''PDF'')(iTr,:), C.(''PDFxi'')(iTr,:)] = DPcalcPDF1D(thisX,cfg.PDF.method,cfg.PDF.xi,1,cfg.PDF.norm,cfg.PDF.kernel,cfg.PDF.width);',...
                    '[C.(''CDF'')(iTr,:), C.(''CDFxi'')(iTr,:)] = DPcalcCDF1D(thisX,cfg.CDF.method,cfg.CDF.xi,1,cfg.CDF.kernel,cfg.CDF.width);',...
                    '[C.(''R'')(iTr,:),C.(''RG'')(iTr,:)] = DPcalcAutoCorr(thisX,cfg.AC.scales,cfg.N,1,cfg.AC.Nsc,cfg.AC.norm,cfg.AC.slopeInds);', ...
                    };

elseif (method==2) %trialTime 
 
     funCommands  = {'[C.(''NE'')(iT,1,iTr), estimation, nbias,sigma] = DPcalcNE(thisX,cfg.NE.Nbins,cfg.NE.approach,cfg.NE.base,1,cfg.NE.norm);', ...
                     'C.(''SE'')(iT,1,iTr) = DPcalcSampEn(thisX,cfg.SE.m,cfg.SE.r,N,1);', ...
                     '[C.(''MSE'')(iT,:,iTr), C.(''STD'')(iT,:,iTr),C.(''G'')(iT,:,iTr)] = DPcalcMSE(thisX,cfg.MSE.m,cfg.MSE.r,cfg.MSE.scales,N,1,cfg.MSE.Nsc,,cfg.MSE.downsample,cfg.MSE.normScale,cfg.MSE.slopeInds);', ...
                     '[C.(''MSEc'')(iT,:,iTr),C.(''RMSEc'')(iT,:,iTr), C.(''STDc'')(iT,:,iTr), C.(''MSEr'')(iT,:,iTr),C.(''RMSEr'')(iT,:,iTr), C.(''STDr'')(iT,:,iTr),C.(''Gc'')(iT,:,iTr),C.(''RGc'')(iT,:,iTr),,C.(''Gr'')(iT,:,iTr),C.(''RGr'')(iT,:,iTr)] = DPcalcMSEall(thisX,cfg.MSEall.m,cfg.MSEall.r,cfg.MSEall.scales,N,1,cfg.MSEall.Nsc,cfg.MSEall.slopeInds);', ...
                     '[C.(''CMSE'')(iT,:,iTr), C.(''CSTD'')(iT,:,iTr),C.(''CG'')(iT,:,iTr)] = DPcalcCMSE(thisX,cfg.CMSE.m,cfg.CMSE.r,cfg.CMSE.scales,N,1,cfg.CMSE.Nsc,cfg.CMSE.normScale,cfg.CMSE.slopeInds);', ...                    
                     '[C.(''H'')(iT,:,iTr), C.(''logF'')(iT,:,iTr) ] = DPcalcDFA(thisX,cfg.DFA.scales,cfg.DFA.logScales,cfg.DFA.order,N,cfg.DFA.Nsc,cfg.DFA.slopeInds);',...
                     '[C.(''P'')(iT,:,iTr), C.(''DOF'')(iT,1,iTr), C.(''logP'')(iT,:,iTr), C.(''PS'')(iT,:,iTr)] = DPcalcSP(thisX,cfg.SP.logf,cfg.SP.winfun,cfg.SP.NFFT,cfg.SP.Nf,N,cfg.SP.slopeInds);',...                     
                     '[C.(''logV'')(iT,:,iTr), C.(''VS'')(iT,:,iTr) ] = DPcalcVarGram(thisX,cfg.VGR.scales,cfg.VGR.logScales,cfg.VGR.Nsc,cfg.VGR.slopeInds);',...
                     '[C.(''PDF'')(iT,:,iTr), C.(''PDFxi'')(iT,:,iTr)] = DPcalcPDF1D(thisX,cfg.PDF.method,cfg.PDF.xi,1,cfg.PDF.norm,cfg.PDF.kernel,cfg.PDF.width);',...
                     '[C.(''CDF'')(iT,:,iTr), C.(''CDFxi'')(iT,:,iTr)] = DPcalcCDF1D(thisX,cfg.CDF.method,cfg.CDF.xi,1,cfg.CDF.kernel,cfg.CDF.width);',...
                     '[C.(''R'')(iT,:,iTr), C.(''RG'')(iT,:,iTr)] = DPcalcAutoCorr(thisX,cfg.AC.scales,N,1,cfg.AC.Nsc,cfg.AC.norm,cfg.AC.slopeInds);', ...
                     };

                
else %ensemble 
    
    funCommands  = {'[C.(''NE'')(iT), estimation, nbias,sigma] = DPcalcNEens(thisX,cfg.NE.Nbins,cfg.NE.approach,cfg.NE.base,cfg.NE.norm);', ...
                    'C.(''SE'')(iT) = DPcalcSampEnEns(thisX,cfg.SE.m,cfg.SE.r,N,cfg.Ntr);', ...
                    '[C.(''MSE'')(iT,:), C.(''STD'')(iT,:),C.(''G'')(iT,:)] = DPcalcMSEens(thisX,cfg.MSE.m,cfg.MSE.r,cfg.MSE.scales,N,cfg.Ntr,cfg.MSE.Nsc,cfg.MSE.downsample,cfg.MSE.normScale,cfg.MSE.slopeInds);', ...
                    '[C.(''MSEc'')(iT,:),C.(''RMSEc'')(iT,:), C.(''STDc'')(iT,:),C.(''MSEr'')(iT,:),C.(''RMSEr'')(iT,:), C.(''STDr'')(iT,:),C.(''Gc'')(iT,:),C.(''RGc'')(iT,:),C.(''Gr'')(iT,:),C.(''RGr'')(iT,:)] = DPcalcMSEallEns(thisX,cfg.MSEall.m,cfg.MSEall.r,cfg.MSEall.scales,N,cfg.Ntr,cfg.MSEall.Nsc,cfg.MSEall.slopeInds);', ...
                    '[C.(''CMSE'')(iT,:), C.(''CSTD'')(iT,:),C.(''CG'')(iT,:)] = DPcalcCMSEens(thisX,cfg.CMSE.m,cfg.CMSE.r,cfg.CMSE.scales,N,cfg.Ntr,cfg.CMSE.Nsc,cfg.CMSE.normScale,cfg.CMSE.slopeInds);', ...
                    '[C.(''H'')(iT,:), C.(''logF'')(iT,:) ] = DPcalcDFAens(thisX,cfg.DFA.scales,cfg.DFA.logScales,cfg.DFA.order,N,cfg.Ntr,cfg.DFA.Nsc,cfg.DFA.slopeInds);',...
                    '[C.(''P'')(iT,:), C.(''DOF'')(iT), C.(''logP'')(iT,:), C.(''PS'')(iT,:)] = DPcalcSPens(thisX,cfg.SP.logf,cfg.SP.winfun,cfg.SP.NFFT,cfg.SP.Nf,N,cfg.SP.slopeInds);',...  
                    '[C.(''logV'')(iT,:), C.(''VS'')(iT,:)] = DPcalcVarGramEns(thisX,cfg.VGR.scales,cfg.VGR.logScales,cfg.VGR.Nsc,cfg.VGR.slopeInds);',...  
                    '[C.(''PDF'')(iT,:), C.(''PDFxi'')(iT,:)] = DPcalcPDF1Dens(thisX,cfg.PDF.method,cfg.PDF.xi,cfg.PDF.norm,cfg.PDF.kernel,cfg.PDF.width);',...
                    '[C.(''CDF'')(iT,:), C.(''CDFxi'')(iT,:)] = DPcalcCDF1Dens(thisX,cfg.CDF.method,cfg.CDF.xi,cfg.CDF.kernel,cfg.CDF.width);',...
                    '[C.(''R'')(iT,:), C.(''RG'')(iT,:)] = DPcalcAutoCorrEns(thisX,cfg.AC.scales,N,cfg.Ntr,cfg.AC.Nsc,cfg.AC.norm,cfg.AC.slopeInds);', ...
                     };
end   


%Number of data points per calculation:
if method==1  %trial
    NpointsPcalc = cfg.N;
else
    NpointsPcalc = cfg.Nwin;  
end
 
%Measures' parameters structures:
if measures(1)
    
    if method==3  %ensemble
        NpointsPcalcNE = NpointsPcalc*cfg.Ntr;
    else
        NpointsPcalcNE = NpointsPcalc;
    end

    %Default number of bins:
    defaultNbins = max(round(sqrt(NpointsPcalcNE)),3);
    
    
    if isfield(cfg,'NE')
        varName='NE';
        testNE = {@(NE)isstruct(NE) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.NE] = DPvalidateData(cfg.NE,testNE,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.NE,'Nbins')
            
            varName='Nbins';
            testNBINS = {@(Nbins)isnumeric(Nbins),...
                         @(Nbins)isreal(Nbins),...
                         @(Nbins)isscalar(Nbins),...
                         @(Nbins,NpointsPcalcNE)( Nbins>2) && ( Nbins<NpointsPcalcNE ),...
                         @(Nbins)round(Nbins)==Nbins };
            param={{},{},{},{NpointsPcalcNE},{}};
            mode=['e','e','e','e','w'];
            execfun={{},{},{},@(Nbins)Nbins(1),@(Nbins)round(Nbins)};
            default=defaultNbins;
            [RESULT, cfg.NE.Nbins] = DPvalidateData(cfg.NE.Nbins,testNBINS,param,mode,execfun,default,varName,funcName);
        else
            cfg.NE.Nbins = defaultNbins;
            %cprintf('Magenta','WARNING: Structure field cfg.NE.Nbins is missing. Replacing with default value: cfg.NE.Nbins = %d.\n',cfg.Nbins)
            fprintf('WARNING: Structure field cfg.NE.Nbins is missing. Replacing with default value: cfg.NE.Nbins = %d.\n',cfg.Nbins)
        end
        
        if isfield(cfg.NE,'approach')
            varName='approach';
            testAPPROACH = { @(approach)ischar(approach),...
                             @(approach)isvector(approach),...
                             @(approach)any(strcmpi(approach,{'unbiased','mmse','biased'})) };
            param={{},{},{}};
            mode=['e','e','e']; 
            execfun={{},{},{}};
            default='biased';
            [RESULT, cfg.NE.approach] = DPvalidateData(cfg.NE.approach,testAPPROACH,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.NE.approach is missing. Replacing with default value cfg.NE.approach=''biased''.\n')
            fprintf('WARNING: Structure field cfg.NE.approach is missing. Replacing with default value cfg.NE.approach=''biased''.\n')
            cfg.NE.approach='biased';
        end

        if isfield(cfg.NE,'base')
            varName='NE.base';
            testBASE = {@(base)isnumeric(base),...
                        @(base)isreal(base),...
                        @(base)isscalar(base),...
                        @(base)(base>0) };
            param={{},{},{},{}};
            mode=['e','e','w','e'];
            execfun={{},{},@(base)base(1),{}};
            default=exp(1);
            [RESULT, cfg.NE.base] = DPvalidateData(cfg.NE.base,testBASE,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.NE.base is missing. Replacing with default value cfg.NE.base=e.\n')
            fprintf('WARNING: Structure field cfg.NE.base is missing. Replacing with default value cfg.NE.base=e.\n')
            cfg.NE.base=exp(1);
        end
        
         if isfield(cfg.NE,'norm')
            varName='norm';
            testNORM = { @(norm)any(norm==[0 1]) };
            param={{}};
            mode=['e']; 
            execfun={{}};
            default='biased';
            [RESULT, cfg.NE.norm] = DPvalidateData(cfg.NE.norm,testNORM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.NE.norm is missing. Replacing with default value cfg.NE.norm=1.\n')
            fprintf('WARNING: Structure field cfg.NE.norm is missing. Replacing with default value cfg.NE.norm=1.\n')
            cfg.NE.norm=1;
         end
        
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.NE is missing. Replacing with default values cfg.NE.Nbins = sqrt(N), cfg.NE.approach=''biased'', cfg.NE.base=e, cfg.NE.norm=1.\n')
        fprintf('WARNING: Structure field cfg.NE is missing. Replacing with default values cfg.NE.Nbins = sqrt(N), cfg.NE.approach=''biased'', cfg.NE.base=e, cfg.NE.norm=1.\n')
        cfg.NE.approach='biased';
        cfg.NE.base=exp(1);
        cfg.NE.Nbins=defaultNbins;
        cfg.NE.norm=1;
    end
    
    %Set the norm flag:
    if (method==1)
        if strcmpi(cfg.normal,'zscore')
            cfg.NE.normInp=0;%...if we haven already normalized trial data...
        else
            cfg.NE.normInp=1;%...if we haven't already normalized trial data
        end
    end
    
end

if measures(2)
    if isfield(cfg,'SE')
        varName='SE';
        testSE = {@(SE)isstruct(SE) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.SE] = DPvalidateData(cfg.SE,testSE,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.SE,'m')
            varName='SE.m';
            testM = {@(m)isnumeric(m),...
                     @(m)isreal(m),...
                     @(m)isscalar(m),...
                     @(m,N)( (m>0)&&(m<N) ),...
                     @(m)(round(m)==m) };
            param={{},{},{},{cfg.N},{}};
            mode=['e','e','w','e','w'];
            execfun={{},{},@(m)m(1),{},@(m)round(m)};
            default=2;
            [RESULT, cfg.SE.m] = DPvalidateData(cfg.SE.m,testM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.SE.m is missing. Replacing with default value cfg.SE.m=2.\n')
            fprintf('WARNING: Structure field cfg.SE.m is missing. Replacing with default value cfg.SE.m=2.\n')            
            cfg.SE.m=2;
        end
        
        if isfield(cfg.SE,'r')
            varName='SE.r';
            testR = {@(r)isnumeric(r),...
                     @(r)isreal(r),...
                     @(r)isscalar(r),...
                     @(r)( (r>0) && (r<1) ) };
            param={{},{},{},{}};
            mode=['e','e','w','e'];
            execfun={{},{},@(r)r(1),{}};
            default=0.5;
            [RESULT, cfg.SE.r] = DPvalidateData(cfg.SE.r,testR,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.SE.r is missing. Replacing with default value cfg.SE.r=0.5.\n')
            fprintf('WARNING: Structure field cfg.SE.r is missing. Replacing with default value cfg.SE.r=0.5.\n')
            cfg.SE.r=0.5;
        end
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.SE is missing. Replacing with default values cfg.SE.m=2, cfg.SE.r=0.5.\n')
        fprintf('WARNING: Structure field cfg.SE is missing. Replacing with default values cfg.SE.m=2, cfg.SE.r=0.5.\n')
        cfg.SE.m=2;
        cfg.SE.r=0.5;
    end 
    
    %Set the norm flag:
    if (method==1)
        if strcmpi(cfg.normal,'zscore')
            cfg.SE.norm=0;%...if we haven already normalized trial data...
        else
            cfg.SE.norm=1;%...if we haven't already normalized trial data
        end
    end
end


if measures(3)
    
     defaultScales = [1:floor(NpointsPcalc/50)].';
            
    if isfield(cfg,'MSE')
        varName='MSE';
        testMSE = {@(MSE)isstruct(MSE) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.MSE] = DPvalidateData(cfg.MSE,testMSE,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.MSE,'m')
            varName='MSE.m';
            testM = {@(m)isnumeric(m),...
                     @(m)isreal(m),...
                     @(m)isscalar(m),...
                     @(m,N)( (m>0)&&(m<N) ),...
                     @(m)(round(m)==m) };
            param={{},{},{},{cfg.N},{}};
            mode=['e','e','w','e','w'];
            execfun={{},{},@(m)m(1),{},@(m)round(m)};
            default=2;
            [RESULT, cfg.MSE.m] = DPvalidateData(cfg.MSE.m,testM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSE.m is missing. Replacing with default value cfg.MSE.m=2.\n')
            fprintf('WARNING: Structure field cfg.MSE.m is missing. Replacing with default value cfg.MSE.m=2.\n')
            cfg.MSE.m=2;
        end
        
        if isfield(cfg.MSE,'r')
            varName='MSE.r';
            testR = {@(r)isnumeric(r),...
                     @(r)isreal(r),...
                     @(r)isscalar(r),...
                     @(r)( (r>0) && (r<1) ) };
            param={{},{},{},{}};
            mode=['e','e','w','e'];
            execfun={{},{},@(r)r(1),{}};
            default=0.5;
            [RESULT, cfg.MSE.r] = DPvalidateData(cfg.MSE.r,testR,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSE.r is missing. Replacing with default value cfg.MSE.r=0.5.\n')
            fprintf('WARNING: Structure field cfg.MSE.r is missing. Replacing with default value cfg.MSE.r=0.5.\n')
            cfg.MSE.r=0.5;
        end
        
        if isfield(cfg.MSE,'scales')
            varName='MSE.scales';
            testSCALES = {@(scales)isnumeric(scales),...
                          @(scales)all(isreal(scales)),...
                          @(scales)~isscalar(scales),...
                          @(scales,NpointsPcalc)all( (scales>0) & (scales<NpointsPcalc/2) ),...
                          @(scales)all(round(scales)==scales),...
                          @(scales)all(unique(scales)==scales)};
            param={{},{},{},{NpointsPcalc},{},{}};
            mode=['e','e','w','e','w','w'];
            execfun={{},{},@(scales)1:scales,{},@(scales)round(scales),@(scales)unique(scales)};
            default=defaultScales;
            [RESULT, cfg.MSE.scales] = DPvalidateData(cfg.MSE.scales,testSCALES,param,mode,execfun,default,varName,funcName);
            cfg.MSE.scales = cfg.MSE.scales(:);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSE.scales is missing. Replacing with default value cfg.MSE.scales=1:floor(N/50).\n')
            fprintf('WARNING: Structure field cfg.MSE.scales is missing. Replacing with default value cfg.MSE.scales=1:floor(N/50).\n')
            cfg.MSE.scales=defaultScales;
        end
        cfg.MSE.Nsc = length(cfg.MSE.scales);

         if isfield(cfg.MSE,'downsample')
            varName='MSE.downsample';
            testDOWNSAMPLE = {@(downsample)any(downsample==[0,1])};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=0;
            [RESULT, cfg.MSE.downsample] = DPvalidateData(cfg.MSE.downsample,testDOWNSAMPLE,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSE.downsample is missing. Replacing with default value cfg.MSE.downsample=0.\n')
            fprintf('WARNING: Structure field cfg.MSE.downsample is missing. Replacing with default value cfg.MSE.downsample=0.\n')
            cfg.MSE.downsample=0;
         end
        
        if isfield(cfg.MSE,'normScale')
            varName='MSE.normScale';
            testNORMSCALE = {@(normScale)any(normScale==[0,1])};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=0;
            [RESULT, cfg.MSE.normScale] = DPvalidateData(cfg.MSE.normScale,testNORMSCALE,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSE.normScale is missing. Replacing with default value cfg.MSE.normScale=0.\n')
            fprintf('WARNING: Structure field cfg.MSE.normScale is missing. Replacing with default value cfg.MSE.normScale=0.\n')
            cfg.MSE.normScale=0;
        end
        
        if isfield(cfg.MSE,'slopeInds')
            if ~isempty(cfg.MSE.slopeInds)
                varName='MSE.slopeInds';
                testSLOPEINDS = {@(slopeInds)isnumeric(slopeInds),...
                                 @(slopeInds)isreal(slopeInds),...
                                 @(slopeInds)all(slopeInds>0),...
                                 @(slopeInds)isvector(slopeInds),...
                                 @(slopeInds)all(round(slopeInds)==slopeInds),...
                                 @(slopeInds)all(unique(slopeInds)==slopeInds)};
                param={{},{},{},{},{},{}};
                mode=['e','e','e','e','e','e'];
                execfun={{},{},{},{},{},{}};
                default=1:cfg.MSE.Nsc;
                [RESULT, cfg.MSE.slopeInds] = DPvalidateData(cfg.MSE.slopeInds,testSLOPEINDS,param,mode,execfun,default,varName,funcName);
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSE.slopeInds is missing. Replacing with default value cfg.MSE.slopeInds=1:cfg.MSE.Nsc.\n')
            fprintf('WARNING: Structure field cfg.MSE.slopeInds is missing. Replacing with default value cfg.MSE.slopeInds=1:cfg.MSE.Nsc.\n')
            cfg.MSE.slopeInds=1:cfg.MSE.Nsc;
        end 
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.MSE is missing. Replacing with default values cfg.MSE.m=2, cfg.MSE.r=0.5, cfg.MSE.scales=[1:floor(N/50)], cfg.MSE.normScale=0.\n')
        fprintf('WARNING: Structure field cfg.MSE is missing. Replacing with default values cfg.MSE.m=2, cfg.MSE.r=0.5, cfg.MSE.scales=[1:floor(N/50)], cfg.MSE.normScale=0.\n')
        cfg.MSE.m=2;
        cfg.MSE.r=0.5;
        cfg.MSE.scales = defaultScales; 
        cfg.MSE.Nsc = length(cfg.MSE.scales);
        cfg.MSE.downsample=0;
        cfg.MSE.normScale=0;
        cfg.MSE.slopeInds=1:cfg.MSE.Nsc;
    end

    if cfg.MSE.Nsc<2
        %cprintf('Magenta','WARNING: Not enough data for MSE. MSE is excluded from calculation.\n')
        fprintf('WARNING: Not enough data for MSE. MSE is excluded from calculation.\n')
        measures(3)=0;
    end
      
    %Set the number of outpus equal to the number of scales for MSE and STD
    %                    MSE         STD      G
    NoutsPmeas{3} = [cfg.MSE.Nsc cfg.MSE.Nsc, 1];
    if (cfg.MSE.Nsc>NoutsMAX)
        NoutsMAX = cfg.MSE.Nsc;
    end
end


if measures(4)
    
     defaultScales = [1:floor(NpointsPcalc/50)].';
            
    if isfield(cfg,'MSEall')
        varName='MSEall';
        testMSEALL = {@(MSEall)isstruct(MSEall) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.MSEall] = DPvalidateData(cfg.MSEall,testMSEALL,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.MSEall,'m')
            varName='MSEall.m';
            testM = {@(m)isnumeric(m),...
                     @(m)isreal(m),...
                     @(m)isscalar(m),...
                     @(m,N)( (m>0)&&(m<N) ),...
                     @(m)(round(m)==m) };
            param={{},{},{},{cfg.N},{}};
            mode=['e','e','w','e','w'];
            execfun={{},{},@(m)m(1),{},@(m)round(m)};
            default=2;
            [RESULT, cfg.MSEall.m] = DPvalidateData(cfg.MSEall.m,testM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSEall.m is missing. Replacing with default value cfg.MSEall.m=2.\n')
            fprintf('WARNING: Structure field cfg.MSEall.m is missing. Replacing with default value cfg.MSEall.m=2.\n')
            cfg.MSEall.m=2;
        end
        
        if isfield(cfg.MSEall,'r')
            varName='MSEall.r';
            testR = {@(r)isnumeric(r),...
                     @(r)isreal(r),...
                     @(r)isscalar(r),...
                     @(r)( (r>0) && (r<1) ) };
            param={{},{},{},{}};
            mode=['e','e','w','e'];
            execfun={{},{},@(r)r(1),{}};
            default=0.5;
            [RESULT, cfg.MSEall.r] = DPvalidateData(cfg.MSEall.r,testR,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSEall.r is missing. Replacing with default value cfg.MSEall.r=0.5.\n')
            fprintf('WARNING: Structure field cfg.MSEall.r is missing. Replacing with default value cfg.MSEall.r=0.5.\n')
            cfg.MSEall.r=0.5;
        end
        
        if isfield(cfg.MSEall,'scales')
            varName='MSEall.scales';
            testSCALES = {@(scales)isnumeric(scales),...
                          @(scales)all(isreal(scales)),...
                          @(scales)~isscalar(scales),...
                          @(scales,NpointsPcalc)all( (scales>0) & (scales<NpointsPcalc/2) ),...
                          @(scales)all(round(scales)==scales),...
                          @(scales)all(unique(scales)==scales)};
            param={{},{},{},{NpointsPcalc},{},{}};
            mode=['e','e','w','e','w','w'];
            execfun={{},{},@(scales)1:scales,{},@(scales)round(scales),@(scales)unique(scales)};
            default=defaultScales;
            [RESULT, cfg.MSEall.scales] = DPvalidateData(cfg.MSEall.scales,testSCALES,param,mode,execfun,default,varName,funcName);
            cfg.MSEall.scales = cfg.MSEall.scales(:);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSEall.scales is missing. Replacing with default value cfg.MSEall.scales=1:floor(N/50).\n')
            fprintf('WARNING: Structure field cfg.MSEall.scales is missing. Replacing with default value cfg.MSEall.scales=1:floor(N/50).\n')
            cfg.MSEall.scales=defaultScales;
        end
        cfg.MSEall.Nsc = length(cfg.MSEall.scales);

        if isfield(cfg.MSEall,'slopeInds')
            if ~isempty(cfg.MSEall.slopeInds)
                varName='MSEall.slopeInds';
                testSLOPEINDS = {@(slopeInds)isnumeric(slopeInds),...
                                 @(slopeInds)isreal(slopeInds),...
                                 @(slopeInds)all(slopeInds>0),...
                                 @(slopeInds)isvector(slopeInds),...
                                 @(slopeInds)all(round(slopeInds)==slopeInds),...
                                 @(slopeInds)all(unique(slopeInds)==slopeInds)};
                param={{},{},{},{},{},{}};
                mode=['e','e','e','e','e','e'];
                execfun={{},{},{},{},{},{}};
                default=1:cfg.MSEall.Nsc;
                [RESULT, cfg.MSEall.slopeInds] = DPvalidateData(cfg.MSEall.slopeInds,testSLOPEINDS,param,mode,execfun,default,varName,funcName);
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSEall.slopeInds is missing. Replacing with default value cfg.MSEall.slopeInds=1:cfg.MSEall.Nsc.\n')
            fprintf('WARNING: Structure field cfg.MSEall.slopeInds is missing. Replacing with default value cfg.MSEall.slopeInds=1:cfg.MSEall.Nsc.\n')
            cfg.MSEall.slopeInds=1:cfg.MSEall.Nsc;
        end
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.MSE is missing. Replacing with default values cfg.MSE.m=2, cfg.MSE.r=0.5, cfg.MSE.scales=[1:floor(N/50)].\n')
        fprintf('WARNING: Structure field cfg.MSEall is missing. Replacing with default values cfg.MSEall.m=2, cfg.MSEall.r=0.5, cfg.MSEall.scales=[1:floor(N/50)].\n')
        cfg.MSEall.m=2;
        cfg.MSEall.r=0.5;
        cfg.MSEall.scales = defaultScales;
        cfg.MSEall.Nsc = length(cfg.MSEall.scales);
        cfg.MSEall.normScale=0;
        cfg.MSEall.slopeInds=1:cfg.MSE.Nsc;
    end

    if cfg.MSEall.Nsc<2
        %cprintf('Magenta','WARNING: Not enough data for MSEall. MSEall is excluded from calculation.\n')
        fprintf('WARNING: Not enough data for MSEall. MSEall is excluded from calculation.\n')
        measures(4)=0;
    end
      
    %Set the number of outputs equal to the number of scales for MSEall and STD
    NoutsPmeas{4} = [cfg.MSEall.Nsc*ones(1,6),ones(1,4)];
    if (cfg.MSEall.Nsc>NoutsMAX)
        NoutsMAX = cfg.MSEall.Nsc;
    end
end


if measures(5)
    
    defaultScales = [1:floor(NpointsPcalc/50)].';
    
    if isfield(cfg,'CMSE')
        varName='CMSE';
        testCMSE = {@(CMSE)isstruct(CMSE) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.CMSE] = DPvalidateData(cfg.CMSE,testCMSE,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.CMSE,'m')
            varName='CMSE.m';
            testM = {@(m)isnumeric(m),...
                     @(m)isreal(m),...
                     @(m)isscalar(m),...
                     @(m,N)( (m>0)&&(m<N) ),...
                     @(m)(round(m)==m) };
            param={{},{},{},{cfg.N},{}};
            mode=['e','e','w','e','w'];
            execfun={{},{},@(m)m(1),{},@(m)round(m)};
            default=2;
            [RESULT, cfg.CMSE.m] = DPvalidateData(cfg.CMSE.m,testM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.CMSE.m is missing. Replacing with default value cfg.CMSE.m=2.\n')
            fprintf('WARNING: Structure field cfg.CMSE.m is missing. Replacing with default value cfg.CMSE.m=2.\n')
            cfg.CMSE.m=2;
        end
        
        if isfield(cfg.CMSE,'r')
            varName='CMSE.r';
            testR = {@(r)isnumeric(r),...
                     @(r)isreal(r),...
                     @(r)isscalar(r),...
                     @(r)( (r>0) && (r<1) ) };
            param={{},{},{},{}};
            mode=['e','e','w','e'];
            execfun={{},{},@(r)r(1),{}};
            default=0.5;
            [RESULT, cfg.CMSE.r] = DPvalidateData(cfg.CMSE.r,testR,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.CMSE.r is missing. Replacing with default value cfg.CMSE.r=0.5.\n')
            fprintf('WARNING: Structure field cfg.CMSE.r is missing. Replacing with default value cfg.CMSE.r=0.5.\n')
            cfg.CMSE.r=0.5;
        end
        
        if isfield(cfg.CMSE,'scales')
            varName='CMSE.scales';
            testSCALES = {@(scales)isnumeric(scales),...
                          @(scales)all(isreal(scales)),...
                          @(scales)~isscalar(scales),...
                          @(scales,NpointsPcalc)all( (scales>0) & (scales<=NpointsPcalc/50) ),...
                          @(scales)all(round(scales)==scales),...
                          @(scales)all(unique(scales)==scales)};
            param={{},{},{},{NpointsPcalc},{},{}};
            mode=['e','e','w','e','w','w'];
            execfun={{},{},@(scales)1:scales,{},@(scales)round(scales),@(scales)unique(scales)};
            default=defaultScales;
            [RESULT, cfg.CMSE.scales] = DPvalidateData(cfg.CMSE.scales,testSCALES,param,mode,execfun,default,varName,funcName);
            cfg.CMSE.scales = cfg.CMSE.scales(:);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.CMSE.scales is missing. Replacing with default value cfg.CMSE.scales=1:floor(N/50).\n')
            fprintf('WARNING: Structure field cfg.CMSE.scales is missing. Replacing with default value cfg.CMSE.scales=1:floor(N/50).\n')
            cfg.CMSE.scales=defaultScales;
        end
        cfg.CMSE.Nsc = length(cfg.CMSE.scales);
        
        if isfield(cfg.CMSE,'normScale')
            varName='CMSE.normScale';
            testNORMSCALE = {@(normScale)any(normScale==[0,1])};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=0;
            [RESULT, cfg.CMSE.normScale] = DPvalidateData(cfg.CMSE.normScale,testNORMSCALE,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.CMSE.normScale is missing. Replacing with default value cfg.CMSE.normScale=0.\n')
            fprintf('WARNING: Structure field cfg.CMSE.normScale is missing. Replacing with default value cfg.CMSE.normScale=0.\n')
            cfg.CMSE.normScale=0;
        end
        
         if isfield(cfg.CMSE,'slopeInds')
            if ~isempty(cfg.CMSE.slopeInds)
                varName='CMSE.slopeInds';
                testSLOPEINDS = {@(slopeInds)isnumeric(slopeInds),...
                                 @(slopeInds)isreal(slopeInds),...
                                 @(slopeInds)all(slopeInds>0),...
                                 @(slopeInds)isvector(slopeInds),...
                                 @(slopeInds)all(round(slopeInds)==slopeInds),...
                                 @(slopeInds)all(unique(slopeInds)==slopeInds)};
                param={{},{},{},{},{},{}};
                mode=['e','e','e','e','e','e'];
                execfun={{},{},{},{},{},{}};
                default=1:cfg.CMSE.Nsc;
                [RESULT, cfg.CMSE.slopeInds] = DPvalidateData(cfg.CMSE.slopeInds,testSLOPEINDS,param,mode,execfun,default,varName,funcName);
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.MSE.slopeInds is missing. Replacing with default value cfg.CMSE.slopeInds=1:cfg.CMSE.Nsc.\n')
            fprintf('WARNING: Structure field cfg.CMSE.slopeInds is missing. Replacing with default value cfg.CMSE.slopeInds=1:cfg.CMSE.Nsc.\n')
            cfg.CMSE.slopeInds=1:cfg.CMSE.Nsc;
        end
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.CMSE is missing. Replacing with default values cfg.CMSE.m=2, cfg.CMSE.r=0.5, cfg.CMSE.scales=[1:floor(N/50)], cfg.CMSE.normScale=0.\n')
        fprintf('WARNING: Structure field cfg.CMSE is missing. Replacing with default values cfg.CMSE.m=2, cfg.CMSE.r=0.5, cfg.CMSE.scales=[1:floor(N/50)], cfg.CMSE.normScale=0.\n')
        cfg.CMSE.m=2;
        cfg.CMSE.r=0.5;
        cfg.CMSE.scales = defaultScales; 
        cfg.CMSE.Nsc = length(cfg.CMSE.scales);
        cfg.CMSE.normScale=0;
        cfg.CMSE.slopeInds=1:cfg.CMSE.Nsc;
    end
    
    if cfg.CMSE.Nsc<2
        %cprintf('Magenta','WARNING: Not enough data for CMSE. CMSE is excluded from calculation.\n')
        fprintf('WARNING: Not enough data for CMSE. CMSE is excluded from calculation.\n')
        measures(5)=0;
    end
    
    %Set the number of outpus equal to the number of scales for CMSE and CSTD
    %                    CMSE         CSTD      CG
    NoutsPmeas{5} = [cfg.CMSE.Nsc cfg.CMSE.Nsc, 1];
    if (cfg.CMSE.Nsc>NoutsMAX)
        NoutsMAX = cfg.CMSE.Nsc;
    end
end


if measures(6)
    
    if isfield(cfg,'DFA')
        varName='DFA';
        testDFA = {@(DFA)isstruct(DFA) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.DFA] = DPvalidateData(cfg.DFA,testDFA,param,mode,execfun,default,varName,funcName);
        
         if isfield(cfg.DFA,'order')
            varName='DFA.order';
            testORDER = {@(order)isnumeric(order),...
                         @(order)isreal(order),...
                         @(order)isscalar(order),...
                         @(order)(order>0),...
                         @(order)(round(order)==order) };
            param={{},{},{},{},{}};
            mode=['e','e','w','e','w'];
            execfun={{},{},@(order)order(1),{},@(order)round(order)};
            default=2;
            [RESULT, cfg.DFA.order] = DPvalidateData(cfg.DFA.order,testORDER,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.DFA.order is missing. Replacing with default value cfg.DFA.order=2.\n')
            fprintf('WARNING: Structure field cfg.DFA.order is missing. Replacing with default value cfg.DFA.order=2.\n')
            cfg.DFA.order=2;
         end
            
        defaultScales = [cfg.DFA.order+2:floor(NpointsPcalc/10)].';
        if isfield(cfg.DFA,'scales')
            varName='DFA.scales';
            testSCALES = {@(scales)isnumeric(scales),...
                          @(scales)all(isreal(scales)),...
                          @(scales)~isscalar(scales),...
                          @(scales,order,NpointsPcalc)all( (scales>order+1) & (scales<=NpointsPcalc/10) ),...
                          @(scales)all(round(scales)==scales),...
                          @(scales)all(unique(scales)==scales)};
            param={{},{},{},{cfg.DFA.order,NpointsPcalc},{},{}};
            mode=['e','e','w','e','w','w'];
            execfun={{},{},@(scales)1:scales,{},@(scales)round(scales),@(scales)unique(scales)};
            default=defaultScales;
            [RESULT, cfg.DFA.scales] = DPvalidateData(cfg.DFA.scales,testSCALES,param,mode,execfun,default,varName,funcName);
            cfg.DFA.scales = cfg.DFA.scales(:);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.DFA.scales is missing. Replacing with default value cfg.DFA.scales=cfg.DFA.order+2:floor(N/10).\n')
            fprintf('WARNING: Structure field cfg.DFA.scales is missing. Replacing with default value cfg.DFA.scales=cfg.DFA.order+2:floor(N/10).\n')
            cfg.DFA.scales=defaultScales;
        end
        cfg.DFA.Nsc = length(cfg.DFA.scales);
        
        if isfield(cfg.DFA,'slopeInds')
            if ~isempty(cfg.DFA.slopeInds)
                varName='DFA.slopeInds';
                testSLOPEINDS = {@(slopeInds)isnumeric(slopeInds),...
                                 @(slopeInds)isreal(slopeInds),...
                                 @(slopeInds)all(slopeInds>0),...
                                 @(slopeInds)isvector(slopeInds),...
                                 @(slopeInds)all(round(slopeInds)==slopeInds),...
                                 @(slopeInds)all(unique(slopeInds)==slopeInds)};
                param={{},{},{},{},{},{}};
                mode=['e','e','e','e','e','e'];
                execfun={{},{},{},{},{},{}};
                default=1:cfg.DFA.Nsc;
                [RESULT, cfg.DFA.slopeInds] = DPvalidateData(cfg.DFA.slopeInds,testSLOPEINDS,param,mode,execfun,default,varName,funcName);
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.DFA.slopeInds is missing. Replacing with default value cfg.DFA.slopeInds=1:cfg.DFA.Nsc.\n')
            fprintf('WARNING: Structure field cfg.DFA.slopeInds is missing. Replacing with default value cfg.DFA.slopeInds=1:cfg.DFA.Nsc.\n')
            cfg.DFA.slopeInds=1:cfg.DFA.Nsc;
        end 
    else
        %cprintf('Magenta','WARNING: Structure field cfg.DFA is missing. Replacing with default values cfg.DFA.order=2, cfg.DFA.scales=[cfg.DFA.order+2:floor(N/10)].\n')
        fprintf('WARNING: Structure field cfg.DFA is missing. Replacing with default values cfg.DFA.order=2, cfg.DFA.scales=[cfg.DFA.order+2:floor(N/10)].\n')
        cfg.DFA.order=2;
        cfg.DFA.scales = [cfg.DFA.order+2:floor(NpointsPcalc/10)].';
        cfg.DFA.Nsc = length(cfg.DFA.scales);
        cfg.DFA.slopeInds=1:cfg.DFA.Nsc;
    end
    cfg.DFA.logScales = log(cfg.DFA.scales*cfg.Ts);

    if cfg.DFA.Nsc<2
        %cprintf('Magenta','WARNING: Not enough data for DFA. DFA is excluded from calculation.\n')
        fprintf('WARNING: Not enough data for DFA. DFA is excluded from calculation.\n')
        measures(6)=0;
    end
    
    %Set the number of outpus equal to the number of scales for DFA
    %                H     logF
    NoutsPmeas{6} = [1 cfg.DFA.Nsc];
    if (cfg.DFA.Nsc>NoutsMAX)
        NoutsMAX = cfg.DFA.Nsc;
    end

end


if measures(7)
   
    %Default NFFT:
    defaultNFFT = 2^nextpow2(NpointsPcalc);
    if isfield(cfg,'SP')
        
        if isfield(cfg.SP,'NFFT')
            varName='SP.NFFT';
            testNFFT = {@(NFFT)isnumeric(NFFT),...
                        @(NFFT)isreal(NFFT),...
                        @(NFFT)isscalar(NFFT),...
                        @(NFFT,defaultNFFT)( (NFFT>2)&&(NFFT<=defaultNFFT) ),...
                        @(NFFT)(round(log2(NFFT))==log2(NFFT)) };
            param={{},{},{},{defaultNFFT},{}};
            mode=['e','e','w','e','e'];
            execfun={{},{},@(NFFT)NFFT(1),{},{}};
            default=defaultNFFT;
            [RESULT, cfg.SP.NFFT] = DPvalidateData(cfg.SP.NFFT,testNFFT,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.SP.NFFT is missing. Replacing with default value cfg.SP.NFFT=2^nextpow2(NpointsPcalc).\n')
            fprintf('WARNING: Structure field cfg.SP.NFFT is missing. Replacing with default value cfg.SP.NFFT=2^nextpow2(NpointsPcalc).\n')
            cfg.SP.NFFT=defaultNFFT;
        end
        cfg.SP.Nf=cfg.SP.NFFT/2;
        
        if isfield(cfg.SP,'winfun')
            varName='SP.winfun';
            testWINFUN = {@(winfun) isa(winfun,'function_handle')};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=@(N)boxcar(N);
            [RESULT, cfg.SP.winfun] = DPvalidateData(cfg.SP.winfun,testWINFUN,param,mode,execfun,default,varName,funcName);
            
        else
           %cprintf('Magenta','WARNING: Structure field cfg.SP.winfun is missing. Applying NO window function.\n')
           fprintf('WARNING: Structure field cfg.SP.win is missing.Replacing with default value cfg.SP.winfun=@(N)boxcar(N).\n')
        end
        
        if isfield(cfg.SP,'slopeInds')
            if ~isempty(cfg.SP.slopeInds)
                varName='SP.slopeInds';
                testSLOPEINDS = {@(slopeInds)isnumeric(slopeInds),...
                                 @(slopeInds)isreal(slopeInds),...
                                 @(slopeInds)all(slopeInds>0),...
                                 @(slopeInds)isvector(slopeInds),...
                                 @(slopeInds)all(round(slopeInds)==slopeInds),...
                                 @(slopeInds)all(unique(slopeInds)==slopeInds)};
                param={{},{},{},{},{},{}};
                mode=['e','e','e','e','e','e'];
                execfun={{},{},{},{},{},{}};
                default=1:cfg.SP.Nf;
                [RESULT, cfg.SP.slopeInds] = DPvalidateData(cfg.SP.slopeInds,testSLOPEINDS,param,mode,execfun,default,varName,funcName);
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.SP.slopeInds is missing. Replacing with default value cfg.SP.slopeInds=1:cfg.SP.Nf.\n')
            fprintf('WARNING: Structure field cfg.SP.slopeInds is missing. Replacing with default value cfg.SP.slopeInds=1:cfg.SP.Nf.\n')
            cfg.SP.slopeInds=1:cfg.SP.Nf;
        end
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.SP is missing. Replacing with default value cfg.SP.NFFT=2^nextpow2(NpointsPcalc) and applying NO window function.\n')
        fprintf('WARNING: Structure field cfg.SP is missing. Replacing with default values cfg.SP.NFFT=2^nextpow2(NpointsPcalc) and cfg.SP.winfun=@(N)boxcar(N).\n')
        cfg.SP.NFFT=defaultNFFT;
        cfg.SP.Nf=cfg.SP.NFFT/2;
        cfg.SP.winfun=@(N)boxcar(N);
        cfg.SP.slopeInds=1:cfg.SP.Nf;
    end

    cfg.SP.f = cfg.fs/2*linspace(0,1,cfg.SP.Nf+1).';
    cfg.SP.logf=log(cfg.SP.f(2:end));
    
   %Set the number of outpus equal to the number of scales for SP
   %                     P       DOF     logP       PS
   NoutsPmeas{7} = [ cfg.SP.Nf,   1,  cfg.SP.Nf,    1];
   if (cfg.SP.Nf>NoutsMAX)
        NoutsMAX = cfg.SP.Nf;
   end
   
end
  

if measures(8)
            
    defaultScales = [1:floor(NpointsPcalc/10)].';
    
    if isfield(cfg,'VGR')
        varName='VGR';
        testVGR = {@(VGR)isstruct(VGR) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.VGR] = DPvalidateData(cfg.VGR,testVGR,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.VGR,'scales')
            varName='VGR.scales';
            testSCALES = {@(scales)isnumeric(scales),...
                          @(scales)all(isreal(scales)),...
                          @(scales)~isscalar(scales),...
                          @(scales,NpointsPcalc)all( (scales>0) & (scales<=NpointsPcalc/10) ),...
                          @(scales)all(round(scales)==scales),...
                          @(scales)all(unique(scales)==scales)};
            param={{},{},{},{NpointsPcalc},{},{}};
            mode=['e','e','w','e','w','w'];
            execfun={{},{},@(scales)1:scales,{},@(scales)round(scales),@(scales)unique(scales)};
            default=defaultScales;
            [RESULT, cfg.VGR.scales] = DPvalidateData(cfg.VGR.scales,testSCALES,param,mode,execfun,default,varName,funcName);
            cfg.VGR.scales= cfg.VGR.scales(:);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.VGR.scales is missing. Replacing with default value cfg.VGR.scales=1:floor(N/10).\n')
            fprintf('WARNING: Structure field cfg.VGR.scales is missing. Replacing with default value cfg.VGR.scales=1:floor(N/10).\n')
            cfg.VGR.scales=defaultScales;
        end
        cfg.VGR.Nsc = length(cfg.VGR.scales);
        
        if isfield(cfg.VGR,'slopeInds')
            if ~isempty(cfg.VGR.slopeInds)
                varName='VGR.slopeInds';
                testSLOPEINDS = {@(slopeInds)isnumeric(slopeInds),...
                                 @(slopeInds)isreal(slopeInds),...
                                 @(slopeInds)all(slopeInds>0),...
                                 @(slopeInds)isvector(slopeInds),...
                                 @(slopeInds)all(round(slopeInds)==slopeInds),...
                                 @(slopeInds)all(unique(slopeInds)==slopeInds)};
                param={{},{},{},{},{},{}};
                mode=['e','e','e','e','e','e'];
                execfun={{},{},{},{},{},{}};
                default=1:cfg.VGR.Nsc;
                [RESULT, cfg.VGR.slopeInds] = DPvalidateData(cfg.VGR.slopeInds,testSLOPEINDS,param,mode,execfun,default,varName,funcName);
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.VGR.slopeInds is missing. Replacing with default value cfg.VGR.slopeInds=1:cfg.VGR.Nsc.\n')
            fprintf('WARNING: Structure field cfg.VGR.slopeInds is missing. Replacing with default value cfg.VGR.slopeInds=1:cfg.VGR.Nsc;.\n')
            cfg.VGR.slopeInds=1:cfg.VGR.Nsc;
        end
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.VGR is missing. Replacing with default values cfg.VGR.scales=[1:floor(N/10)].\n')
        fprintf('WARNING: Structure field cfg.VGR is missing. Replacing with default values cfg.VGR.scales=[1:floor(N/10)].\n')
        cfg.VGR.scales = defaultScales; 
        cfg.VGR.Nsc = length(cfg.VGR.scales);
        cfg.VGR.slopeInds=1:cfg.VGR.Nsc;
    end
    
    cfg.VGR.logScales = log(cfg.VGR.scales*cfg.Ts);

    if cfg.VGR.Nsc<2
        %cprintf('Magenta','WARNING: Not enough data for VGR. VGR is excluded from calculation.\n')
        fprintf('WARNING: Not enough data for VGR. VGR is excluded from calculation.\n')
        measures(8)=0;
    end
    
   %Set the number of outputs equal to the number of scales for VGR
   %                    logV     VS
   NoutsPmeas{8} = [cfg.VGR.Nsc, 1];
   if (cfg.VGR.Nsc>NoutsMAX)
        NoutsMAX = cfg.VGR.Nsc;
   end

    
end


if any(measures([9,10])) %for the two measures of density functions
    minX = min(x(:));
    maxX = max(x(:));
    rangeX = maxX - minX;
    
    if method==3  %ensemble
        NpointsPcalcDF = NpointsPcalc*cfg.Ntr;
    else
        NpointsPcalcDF = NpointsPcalc;
    end
    
    %Default number of bins:
    defaultNbinsDF = max(round(sqrt(NpointsPcalcDF)),3);
    
end


if measures(9)
    
    %Default bin...
    %...centers in the signal's space:
    xbinWidth = rangeX/(defaultNbinsDF-1);
    defaultXIx = [minX:xbinWidth:maxX].';
    %...edges in probability space (quantiles):
    %defaultXIp = [0:1/defaultNbinsDF:1].';
    defaultXIp = [0:1/(defaultNbinsDF+2):1].';
    defaultXIp = defaultXIp(2:end-1,1);
    
    if isfield(cfg,'PDF')
        
        if isfield(cfg.PDF,'method')
            varName='PDF.method';
            testMETHOD = { @(method)ischar(method),...
                           @(method)isvector(method),...
                           @(method)any(strcmpi(method,{'isohist','isodist','kernelisohist','kernelisodist'})) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='isodist';
            [RESULT, cfg.PDF.method] = DPvalidateData(cfg.PDF.method,testMETHOD,param,mode,execfun,default,varName,funcName);
        else
            fprintf('WARNING: Structure field cfg.PDF.method is missing. Replacing with default value cfg.PDF.method=''isodist''.\n')
        end
        
        if strfind(cfg.PDF.method,'isodist')
            defaultXI = defaultXIx;
        else
            defaultXI = defaultXIp;
        end
         
            
        if isfield(cfg.PDF,'bins')
            
            varName='PDF.bins';
            
            testBINS = {@(bins)isnumeric(bins),...
                        @(bins)all(isreal(bins)),...
                        @(bins)isvector(bins) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default=defaultXI;
            [RESULT, cfg.PDF.bins] = DPvalidateData(cfg.PDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
            
            if isequal(cfg.PDF.bins,defaultXI) %if bins has defaulted
                
                cfg.PDF.xi=defaultXI;
            
            elseif isscalar(cfg.PDF.bins) %if bins is the number of bins
                
                testBINS = {@(bins)bins>0,...
                            @(bins)bins>=3,...
                            @(bins)round(bins)==bins };
                param={{},{},{}};
                mode=['e','w','w'];
                execfun={{},@(bins)3,@(bins)round(bins)};
                default=defaultNbinsDF;
                [RESULT, cfg.PDF.bins] = DPvalidateData(cfg.PDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
                
                if strfind(cfg.PDF.method,'isodist')
                    cfg.PDF.xi = [minX:rangeX/(cfg.PDF.bins-1):maxX].';
                else
                    cfg.PDF.xi = [minX:rangeX/cfg.PDF.bins:maxX].';
                    
                end
                cfg.PDF.Nbins = cfg.PDF.bins;
                
            else %bins center of x bins or edges of probability bins (quantiles)
                
                if strfind(cfg.PDF.method,'isodist')
                    testBINS = {@(bins)length(bins)>=3,...
                                @(bins)all(unique(bins)==bins)};
                    param={{},{}};
                    mode=['e','w'];
                    execfun={{},@(bins)unique(bins)};
                    default=defaultXI;
                    [RESULT, cfg.PDF.bins] = DPvalidateData(cfg.PDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
                    cfg.PDF.xi = cfg.PDF.bins(:);
                    cfg.PDF.Nbins = length(cfg.PDF.xi);
                else
                    testBINS = {@(bins)length(bins)>=4,...
                                @(bins)all(unique(bins)==bins)
                                @(bins)all( (bins>=0) & (bins<=1) )};
                    param={{},{},{}};
                    mode=['e','w','e'];
                    execfun={{},@(bins)unique(bins),{}};
                    default=defaultXI;
                    [RESULT, cfg.PDF.bins] = DPvalidateData(cfg.PDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
                    cfg.PDF.xi = cfg.PDF.bins(:);
                    cfg.PDF.Nbins = length(cfg.PDF.xi)-1;
                end
                
                
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.PDF.bins is missing. Setting default value for cfg.PDF.xi.\n')
            fprintf('WARNING: Structure field cfg.PDF.bins is missing. Setting default value for cfg.PDF.xi.\n')
            cfg.PDF.xi=defaultXI;
            cfg.PDF.Nbins=defaultNbinsDF;
        end
        
         if isfield(cfg.PDF,'norm')
            varName='PDF.norm';
            testNORM = {@(norm)any(norm==[0,1])};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=1;
            [RESULT, cfg.PDF.norm] = DPvalidateData(cfg.PDF.norm,testNORM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.PDF.norm is missing. Replacing with default value cfg.PDF.norm=1.\n')
            fprintf('WARNING: Structure field cfg.PDF.norm is missing. Replacing with default value cfg.PDF.norm=1.\n')
            cfg.PDF.norm=1;
         end
        
        if strfind(cfg.PDF.method,'kernel')
            if isfield(cfg.PDF,'kernel')
                testKERNEL = { @(kernel)ischar(kernel) || isa(kernel,'function_handle') };
                param={{}};
                mode=['e'];
                execfun={{}};
                default='normal';
                [RESULT, cfg.PDF.kernel] = DPvalidateData(cfg.PDF.kernel,testKERNEL,param,mode,execfun,default,varName,funcName);
            else
               %cprintf('Magenta','WARNING: Structure field cfg.PDF.kernel is missing. Setting default value ''normal''.\n')
               fprintf('WARNING: Structure field cfg.PDF.kernel is missing. Setting default value ''normal''.\n')
               cfg.PDF.kernel = 'normal';
            end
            
            %Default kernel width equals two times the default bin width in
            %x space: 2*xbinWidth
            if isfield(cfg.PDF,'width')
                testWIDTH = {@(width)isnumeric(width),...
                             @(width)isscalar(width),...
                             @(width)isreal(width),...
                             @(width)width>0,...
                             @(width,rangeX)width<=rangeX };
                param={{},{},{},{},{rangeX}};
                mode=['e','e','e','e','e'];
                execfun={{},{},{},{},{}};
                default=xbinWidth;
                [RESULT, cfg.PDF.width] = DPvalidateData(cfg.PDF.width,testWIDTH,param,mode,execfun,default,varName,funcName);
            else
                %cprintf('Magenta','WARNING: Structure field cfg.PDF.width is missing. Setting default value ''normal''.\n')
               fprintf('WARNING: Structure field cfg.PDF.width is missing. Setting default value 2*rangeX/(defaultNbinsDF-1).\n')
               cfg.PDF.width = xbinWidth;
            end
            
        else
            cfg.PDF.kernel = '';
            cfg.PDF.width = [];
        end
    else
        cfg.PDF.method = 'isodist';
        cfg.PDF.xi = defaultXIx;
        cfg.PDF.Nbins = defaultNbinsDF;
        cfg.PDF.norm=1;
        cfg.PDF.kernel = '';
        cfg.PDF.width = [];
    end
    
    
    %Set the number of outpus equal to the number of bins for PDF
    %                    PDF
    NoutsPmeas{9} = [cfg.PDF.Nbins cfg.PDF.Nbins];
    if (cfg.PDF.Nbins>NoutsMAX)
        NoutsMAX = cfg.PDF.Nbins;
    end
end



if measures(10)

    %Default bin...
    %...edges in the signal's space:
    xbinWidth = rangeX/(defaultNbinsDF-1);
    defaultXIx = [minX:rangeX/defaultNbinsDF:maxX].';
    %...edges in probability space (quantiles):
    %defaultXIp = [0:1/defaultNbinsDF:1].';
    defaultXIp = [0:1/(defaultNbinsDF+2):1].';
    defaultXIp = defaultXIp(2:end-1,1);
    
    if isfield(cfg,'CDF')
        
        if isfield(cfg.CDF,'method')
            varName='CDF.method';
            testMETHOD = { @(method)ischar(method),...
                           @(method)isvector(method),...
                           @(method)any(strcmpi(method,{'quantl','interpol','kernelquantl','kernelinterpol'})) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='quantl';
            [RESULT, cfg.CDF.method] = DPvalidateData(cfg.CDF.method,testMETHOD,param,mode,execfun,default,varName,funcName);
        else
            fprintf('WARNING: Structure field cfg.CDF.method is missing. Replacing with default value cfg.CDF.method=''quantl''.\n')
        end
        
        if strfind(cfg.CDF.method,'interpol')
            defaultXI = defaultXIx;
        else
            defaultXI = defaultXIp;
        end
         
            
        if isfield(cfg.CDF,'bins')
            
            varName='CDF.bins';
            
            testBINS = {@(bins)isnumeric(bins),...
                        @(bins)all(isreal(bins)),...
                        @(bins)isvector(bins) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default=defaultXI;
            [RESULT, cfg.CDF.bins] = DPvalidateData(cfg.CDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
            
            if isequal(cfg.CDF.bins,defaultXI) %if bins has defaulted
                
                cfg.CDF.xi=defaultXI;
            
            elseif isscalar(cfg.CDF.bins) %if bins is the number of bins
                
                testBINS = {@(bins)bins>0,...
                            @(bins)bins>=3,...
                            @(bins)round(bins)==bins };
                param={{},{},{}};
                mode=['e','w','w'];
                execfun={{},@(bins)3,@(bins)round(bins)};
                default=defaultNbinsDF;
                [RESULT, cfg.CDF.bins] = DPvalidateData(cfg.CDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
                
                if strfind(cfg.CDF.method,'interpol')
                    cfg.CDF.xi = [minX:rangeX/(cfg.CDF.bins-1):maxX].';
                else
                    cfg.CDF.xi = [minX:rangeX/cfg.CDF.bins:maxX].';
                    
                end
                cfg.CDF.Nbins = cfg.CDF.bins;
                
            else %bins center of x bins or edges of probability bins (quantiles)
                
                if strfind(cfg.CDF.method,'interpol')
                    testBINS = {@(bins)length(bins)>=4,...
                                @(bins)all(unique(bins)==bins)};
                    param={{},{}};
                    mode=['e','w'];
                    execfun={{},@(bins)unique(bins)};
                    default=defaultXI;
                    [RESULT, cfg.CDF.bins] = DPvalidateData(cfg.CDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
                    cfg.CDF.xi = cfg.CDF.bins(:);
                    cfg.CDF.Nbins = length(cfg.CDF.xi)-1;
                else
                    testBINS = {@(bins)length(bins)>=4,...
                                @(bins)all(unique(bins)==bins)
                                @(bins)all( (bins>=0) & (bins<=1) )};
                    param={{},{},{}};
                    mode=['e','w','e'];
                    execfun={{},@(bins)unique(bins),{}};
                    default=defaultXI;
                    [RESULT, cfg.CDF.bins] = DPvalidateData(cfg.CDF.bins,testBINS,param,mode,execfun,default,varName,funcName);
                    cfg.CDF.xi = cfg.CDF.bins(:);
                    cfg.CDF.Nbins = length(cfg.CDF.xi)-1;
                end
                
                
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.CDF.bins is missing. Setting default value for cfg.CDF.xi.\n')
            fprintf('WARNING: Structure field cfg.CDF.bins is missing. Setting default value for cfg.CDF.xi.\n')
            cfg.CDF.xi=defaultXI;
            cfg.CDF.Nbins=defaultNbinsDF;
        end
        
        
        if strfind(cfg.CDF.method,'kernel')
            if isfield(cfg.CDF,'kernel')
                testKERNEL = { @(kernel)ischar(kernel) || isa(kernel,'function_handle') };
                param={{}};
                mode=['e'];
                execfun={{}};
                default='normal';
                [RESULT, cfg.CDF.kernel] = DPvalidateData(cfg.CDF.kernel,testKERNEL,param,mode,execfun,default,varName,funcName);
            else
               %cprintf('Magenta','WARNING: Structure field cfg.CDF.kernel is missing. Setting default value ''normal''.\n')
               fprintf('WARNING: Structure field cfg.CDF.kernel is missing. Setting default value ''normal''.\n')
               cfg.CDF.kernel = 'normal';
            end
            
            %Default kernel width equals two times the default bin width in
            %x space: 2*xbinWidth
            if isfield(cfg.CDF,'width')
                testWIDTH = {@(width)isnumeric(width),...
                             @(width)isscalar(width),...
                             @(width)isreal(width),...
                             @(width)width>0,...
                             @(width,rangeX)width<=rangeX };
                param={{},{},{},{},{rangeX}};
                mode=['e','e','e','e','e'];
                execfun={{},{},{},{},{}};
                default=xbinWidth;
                [RESULT, cfg.CDF.width] = DPvalidateData(cfg.CDF.width,testWIDTH,param,mode,execfun,default,varName,funcName);
            else
                %cprintf('Magenta','WARNING: Structure field cfg.CDF.width is missing. Setting default value ''normal''.\n')
               fprintf('WARNING: Structure field cfg.CDF.width is missing. Setting default value 2*rangeX/(defaultNbinsDF-1).\n')
               cfg.CDF.width = xbinWidth;
            end
            
        else
            cfg.CDF.kernel = '';
            cfg.CDF.width = [];
        end
    else
        cfg.CDF.method = 'quantl';
        cfg.CDF.xi = defaultXIp;
        cfg.CDF.Nbins = defaultNbinsDF;
        cfg.CDF.kernel = '';
        cfg.CDF.width = [];
    end
    
    
    %Set the number of outpus equal to the number of bins for CDF
    %                    CDF
    CDFNbins1 = cfg.CDF.Nbins+1;
    NoutsPmeas{10} = [CDFNbins1 CDFNbins1];
    if (CDFNbins1>NoutsMAX)
        NoutsMAX = CDFNbins1;
    end
end

if measures(11)
    
    defaultScales = [1:floor(NpointsPcalc/50)].';
            
    if isfield(cfg,'AC')
        varName='AC';
        testAC = {@(AC)isstruct(AC) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.AC] = DPvalidateData(cfg.AC,testAC,param,mode,execfun,default,varName,funcName);
       
        
        if isfield(cfg.AC,'scales')
            varName='AC.scales';
            testSCALES = {@(scales)isnumeric(scales),...
                          @(scales)all(isreal(scales)),...
                          @(scales)~isscalar(scales),...
                          @(scales,NpointsPcalc)all( (scales>0) & (scales<NpointsPcalc/2) ),...
                          @(scales)all(round(scales)==scales),...
                          @(scales)all(unique(scales)==scales)};
            param={{},{},{},{NpointsPcalc},{},{}};
            mode=['e','e','w','e','w','w'];
            execfun={{},{},@(scales)1:scales,{},@(scales)round(scales),@(scales)unique(scales)};
            default=defaultScales;
            [RESULT, cfg.AC.scales] = DPvalidateData(cfg.AC.scales,testSCALES,param,mode,execfun,default,varName,funcName);
            cfg.AC.scales = cfg.AC.scales(:);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.AC.scales is missing. Replacing with default value cfg.AC.scales=1:floor(N/50).\n')
            fprintf('WARNING: Structure field cfg.AC.scales is missing. Replacing with default value cfg.AC.scales=1:floor(N/50).\n')
            cfg.AC.scales=defaultScales;
        end
        cfg.AC.Nsc = length(cfg.AC.scales);
        
        if isfield(cfg.AC,'norm')
            varName='AC.norm';
            testNORM = {@(norm)any(norm==[0,1])};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=1;
            [RESULT, cfg.AC.norm] = DPvalidateData(cfg.AC.norm,testNORM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.AC.norm is missing. Replacing with default value cfg.AC.norm=0.\n')
            fprintf('WARNING: Structure field cfg.AC.norm is missing. Replacing with default value cfg.AC.norm=0.\n')
            cfg.AC.norm=1;
        end
        
        if isfield(cfg.AC,'slopeInds')
            if ~isempty(cfg.AC.slopeInds)
                varName='AC.slopeInds';
                testSLOPEINDS = {@(slopeInds)isnumeric(slopeInds),...
                                 @(slopeInds)isreal(slopeInds),...
                                 @(slopeInds)all(slopeInds>0),...
                                 @(slopeInds)isvector(slopeInds),...
                                 @(slopeInds)all(round(slopeInds)==slopeInds),...
                                 @(slopeInds)all(unique(slopeInds)==slopeInds)};
                param={{},{},{},{},{},{}};
                mode=['e','e','e','e','e','e'];
                execfun={{},{},{},{},{},{}};
                default=1:cfg.AC.Nsc;
                [RESULT, cfg.AC.slopeInds] = DPvalidateData(cfg.AC.slopeInds,testSLOPEINDS,param,mode,execfun,default,varName,funcName);
            end
        else
            %cprintf('Magenta','WARNING: Structure field cfg.AC.slopeInds is missing. Replacing with default value cfg.AC.slopeInds=1:cfg.AC.Nsc.\n')
            fprintf('WARNING: Structure field cfg.AC.slopeInds is missing. Replacing with default value cfg.AC.slopeInds=1:cfg.AC.Nsc.\n')
            cfg.AC.slopeInds=1:cfg.AC.Nsc;
        end 
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.AC is missing. Replacing with default values cfg.AC.scales=[1:floor(N/50)], cfg.AC.norm=1.\n')
        fprintf('WARNING: Structure field cfg.AC is missing. Replacing with default values cfg.AC.scales=[1:floor(N/50)], cfg.AC.norm=1.\n')
        cfg.AC.scales = defaultScales; 
        cfg.AC.Nsc = length(cfg.AC.scales);
        cfg.AC.norm=1;
        cfg.AC.slopeInds=1:cfg.AC.Nsc;
    end

    if cfg.AC.Nsc<2
        %cprintf('Magenta','WARNING: Not enough scales for AC. AC is excluded from calculation.\n')
        fprintf('WARNING: Not enough scales for AC. AC is excluded from calculation.\n')
        measures(11)=0;
    end
      
    %Set the number of outpus equal to the number of scales for R and RG
    %                    R        RG
    NoutsPmeas{11} = [cfg.AC.Nsc  1];
    if (cfg.AC.Nsc>NoutsMAX)
        NoutsMAX = cfg.AC.Nsc;
    end
end


measInds = find(measures); %Indexes of measures to be calculated
Nmeasures = length(measInds); %number of measures to calculate
if (Nmeasures==0)
    error('Exiting because there are no measures to calculate!')
end



%Surrogate statistics structure:
if isfield(cfg,'stats')
    
    varName='stats';
    testSTATS = { @(stats)isstruct(stats)||isempty(stats) };
    param = {{}};
    mode=['e'];
    execfun={{}};
    default=nan;
    [RESULT, cfg.stats] = DPvalidateData(cfg.stats,testSTATS,param,mode,execfun,default,varName,funcName);
    
    
    if ~isempty(cfg.stats)
        
        if isfield(cfg.stats,'method')
            varName='stats.method';
            testSTATMETHOD = { @(statMethod)ischar(statMethod),...
                               @(statMethod)isvector(statMethod),...
                               @(statMethod)strcmpi(statMethod,{'timeshuffling'})||strcmpi(statMethod,{'phaseshuffling'})||strcmpi(statMethod,{'phaserandom'}) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='timeshuffling';
            [RESULT, cfg.stats.method] = DPvalidateData(cfg.stats.method,testSTATMETHOD,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.stats.method is missing. Setting default value ''trialshuffling''.\n')
            fprintf('WARNING: Structure field cfg.stats.method is missing. Setting default value ''trialshuffling''.\n')
            cfg.stats.method = 'timeshuffling';
        end
        
        if isfield(cfg.stats,'alpha')
            varName='stats.alpha';
            testSTATSALPHA = {@(alpha)isnumeric(alpha),...
                              @(alpha)isreal(alpha),...
                              @(alpha)isscalar(alpha) ,...
                              @(alpha)( alpha>0 ) && (alpha<=0.5)};
            param={{},{},{},{}};
            mode=['e','e','e','e'];
            execfun={{},{},{},{}};
            default=0.05;
            [RESULT, cfg.stats.alpha] = DPvalidateData(cfg.stats.alpha,testSTATSALPHA,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.stats.alpha is missing. Setting default value 0.05.\n')
            fprintf('WARNING: Structure field cfg.stats.alpha is missing. Setting default value 0.05.\n')
            cfg.stats.alpha = 0.05;
        end
        
        if isfield(cfg.stats,'tail')
            varName='stats.tail';
            testSTATSTAIL = {@(tail) (tail==1) || (tail==2)};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=1;
            [RESULT, cfg.stats.tail] = DPvalidateData(cfg.stats.tail,testSTATSTAIL,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.stats.tail is missing. Setting default value 1.\n')
            fprintf('WARNING: Structure field cfg.stats.tail is missing. Setting default value 1.\n')
            cfg.stats.tail = 1;
        end
        
        
        if isfield(cfg.stats,'corrMultComp')
            varName='stats.corrMultComp';
            testSTATMULTCOMP = {@(corrMultComp)ischar(corrMultComp),...
                                @(corrMultComp)isvector(corrMultComp),...
                                @(corrMultComp)any( [strcmpi(corrMultComp,{'BONF'}), strcmpi(corrMultComp,{'FDR'}), strcmpi(corrMultComp,{''})]  )  };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default={''};
            [RESULT, cfg.stats.corrMultComp] = DPvalidateData(cfg.stats.corrMultComp,testSTATMULTCOMP,param,mode,execfun,default,varName,funcName);
        else
            cfg.stats.corrMultComp='';
        end
        
        
        defaultNperm = cfg.stats.tail*ceil(1/cfg.stats.alpha);
        if strcmpi(cfg.stats.corrMultComp,'BONF')
            switch method
                case 1
                    defaultNperm = cfg.Ntr*NoutsMAX*defaultNperm;
                case 2
                    defaultNperm =  cfg.Ntr*cfg.Ncalc*NoutsMAX*defaultNperm;
                case 3
                    defaultNperm =  cfg.Ncalc*NoutsMAX*defaultNperm;
            end
        end
        
        if isfield(cfg.stats,'Nperm')
            varName='stats.Nperm';
            testSTATSNPERM = {@(Nperm)isnumeric(Nperm),...
                              @(Nperm)isreal(Nperm),...
                              @(Nperm)isscalar(Nperm) ,...
                              @(Nperm)ceil(Nperm)==Nperm,...
                              @(Nperm,defaultNperm)( Nperm>=defaultNperm ) };
            param={{},{},{},{},{defaultNperm}};
            mode=['e','e','e','w','w'];
            execfun={{},{},{},@(Nperm)ceil(Nperm),@(Nperm)Nperm};
            default=defaultNperm;
            [RESULT, cfg.stats.Nperm] = DPvalidateData(cfg.stats.Nperm,testSTATSNPERM,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.stats.Nperm is missing. Setting default value tail*ceil(1/alpha)=%f, possibly adjusted for Bonferroni correction for multiple comparisons.\n',defaultNperm)
            fprintf('WARNING: Structure field cfg.stats.Nperm is missing. Setting default value tail*ceil(1/alpha)=%f, possibly adjusted for Bonferroni correction for multiple comparisons.\n',defaultNperm)
            cfg.stats.Nperm = defaultNperm;
        end
        
        
        if isfield(cfg.stats,'pointStatMethod')
            varName='stats.pointStatMethod';
            testPOINTSTATMETHOD = {@(pointStatMethod)ischar(pointStatMethod),...
                                   @(pointStatMethod)strcmpi(pointStatMethod,{'t'})||strcmpi(pointStatMethod,{'z'})||strcmpi(pointStatMethod,{''}) };
            param={{},{}};
            mode=['e','e'];
            execfun={{},{}};
            default='';
            [RESULT, cfg.stats.pointStatMethod] = DPvalidateData(cfg.stats.pointStatMethod,testPOINTSTATMETHOD,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.stats.pointStatMethod is missing. Setting default value ''''.\n')
            fprintf('WARNING: Structure field cfg.stats.pointStatMethod is missing. Setting default value ''''.\n')
            cfg.stats.pointStatMethod = '';
        end
        
        
        if isfield(cfg.stats,'multiStatfun')
            
            varName='stats.multiStatfun';
            testMULTISTATSFUN = {@(multiStatfun) isstruct(multiStatfun)||isempty(multiStatfun)};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=[];
            [RESULT, cfg.stats.multiStatfun] = DPvalidateData(cfg.stats.multiStatfun,testMULTISTATSFUN,param,mode,execfun,default,varName,funcName);
            
            if ~isempty(cfg.stats.multiStatfun)
                
                if isfield(cfg.stats.multiStatfun,'fun')
                    varName='stats.multiStatfun.fun';
                    testMULTISTATSFUNCTION = {@(fun) isa(fun,'function_handle')};
                    param={{}};
                    mode=['e'];
                    execfun={{}};
                    default=nan;
                    [RESULT, cfg.stats.multiStatfun.fun] = DPvalidateData(cfg.stats.multiStatfun.fun,testMULTISTATSFUNCTION,param,mode,execfun,default,varName,funcName);
                    
                else
                    error('Structure field cfg.stats.multiStatfun.fun is missing.')
                end
                
                if ~isfield(cfg.stats.multiStatfun,'params')
                    cfg.stats.multiStatfun.params=[];
                end
            end
        else
            cfg.stats.multiStatfun=[];
        end
        
        
    end
    
else
    cfg.stats=[];
end

        
%Select the function commands of the measures to be calculated
thisfunCommands = funCommands(measInds); 
Ncomnds = length(thisfunCommands);

if ~isempty(cfg.stats)
    %Create a matrix of handles to function of surrogate time series' creation:
    switch cfg.stats.method
        case 'timeshuffling'
            cfg.stats.surrfun = @(Data,D,Ntr)DPsurrShufflTime(Data,D,Ntr);
        case 'phaseshuffling'
            cfg.stats.surrfun = @(Data,D,Ntr)DPsurrShufflPhase(Data,D,Ntr);
        case 'phaserandom'
            cfg.stats.surrfun = @(Data,D,Ntr)DPsurrRandPhase(Data,D,Ntr);
    end
end


