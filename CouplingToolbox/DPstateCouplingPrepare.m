function [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,MeasNames, NmeasPmeas, x,y,z] = DPstateCouplingPrepare(cfg,x,y,z)

%This function prepares data and configuration struture for DPstateCoupling
%

%Inputs: 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-cfg: configuration structure with fields
%   -fs: sampling frequency in herz, positive real scalar
%   -fc: central frequency of the signals in herz, 2 (or 3) element vector of
%        positive real values
%   -D: number of time delay dimensions of x, y (and z), a vector of 2 or 3
%       positive integers
%   -tau: a vector of D-1 embedding time delays of x, y (and z) in secs, 
%       a vector of 2 or 3 positive integers 
%   -time: time vector in secs, real valued vector, default=[0:N-1]/fs
%   -Wth: the Theiler window to be used for the nearest neighbors
%         statistics in secs
%--------------------------------------------------------------------------
%   -method: one of 
%            1. 'trial', for estimation per trial
%            2. 'trialTime', for estimation per trial with a sliding time 
%               window
%            3. 'ensemble', for estimation across trials, either pointwise
%            or with a time window
%--------------------------------------------------------------------------
%   -normal: 'zscore' or 'meanCenter', 'linear', string, default='none' 
%           if normal == 'linear', there should also be a vector of 2 real
%           numbers, normVal for the minimum and maximum value, 
%           default, normVal = [-1 1]
%
%           but only if domain='real' 
%--------------------------------------------------------------------------
%Optionally for 'trialTime' or 'ensemble'
%   -timeCalc: a subset of time (in secs), with the time points where the calculation
%       will be performed, default=time
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
%              1.'TE' for transfer entropy (using Gomez-Herrero et al estimator)
%              2.'MI' for mutual information (using Gomez-Herrero et al estimator)
%              3.'PTE' partial transfer entropy (using Gomez-Herrero et al estimator)
%              4.'PMI' partial transinformation (using Gomez-Herrero et al estimator)
%                 cell of strings, string 'all' for all of them, default:
%                 'all'
%              5.'TI' for transinformation (naive estimator)
%              6.'CCD' for conditional coupling divergence
%			   7.'SL' for synchronization likelihood
%              8.'NI' for nonlinear interdependance
%              9.'MP' for mutual prediction

%   -TE: optional structure of inputs related to TE
%        -u: range of time lags to be tested, in secs, vector of real numbers >=0
%            the test is calculated symmetrically for time lages +/- u(i), 
%            default, u=1/fs, i.e. one sampling point
%       -k: number of nearest neighbors
%   -MI: optional structure of inputs related to MI
%       -k: number of nearest neighbors
%       -u: range of time lags to be tested, in secs, vector of real numbers >=0
%           the test is calculated symmetrically for time lages +/- u(i), 
%           default, u=0,
%   -PTE: optional structure of inputs related to PTE
%       -k: number of nearest neighbors
%       -u: range of time lags to be tested, in secs, vector of real numbers >=0
%           the test is calculated symmetrically for time lages +/- u(i), 
%           default, u=1/fs, i.e. one sampling point
%   -PMI: optional structure of inputs related to PMI
%       -k: number of nearest neighbors
%       -u: range of time lags to be tested, in secs, vector of real numbers >=0
%           the test is calculated symmetrically for time lages +/- u(i), 
%           default, u=0,
%   -TI: optional structure of inputs related to TI
%       -r: range of search
%       -u: range of time lags to be tested, in secs, vector of real numbers >=0
%           the test is calculated symmetrically for time lages +/- u(i), 
%           default, u=0,
%   -CCD: optional structure of inputs related to CCD
%       -r: range of search
%       -l: length of check for conditional coupling divergence in secs
%   -SL: optional structure of inputs related to NI
%       -k: number of nearest neighbors
%       -u: range of time lags to be tested, in secs, vector of real numbers >=0
%           the test is calculated symmetrically for time lages +/- u(i), 
%           default, u=0,
%   -NI: optional structure of inputs related to NI
%       -k: number of nearest neighbors
%       -u: range of time lags to be tested, in secs, vector of real numbers >=0
%           the test is calculated symmetrically for time lages +/- u(i), 
%           default, u=0,
%   -MP: optional structure of inputs related to MP
%       -k: number of neighbors
%       -predFun: a handle to a prediction function
%       -prederrFun: a handle to a prediction performance function
%--------------------------------------------------------------------------
%Optionally if statistics are to be made:
%   -stats: a structure with fields:
%       -method:a string with the surrogate test method that is to be applied:
%                 -"trialshuffling": random trial shuffling
%                default: "trialshuffling"
%       -surrfun: a handle to a function of the form
%                 xSurr = surrfun(x,N,D), where
%                   -x: the data series columnwise matrix, 
%                   -N: the number of data points per dimension of the set,
%                    here it is N=cfg.Np
%                   -D: the number of dimensions of the set, here it is D=2
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
%--------------------------------------------------------------------------
% -x: the first signal, time series of N time points, real valued signal: 
%     a matrix N x Ntr
%     
% -y: the second signal, time series of N time points, real valued signal: 
%     a matrix N x Ntr
%
% Optionally for partial measures:
% -z: the third signal for calculation of partial measures,
%     real valued signal: a matrix N x Ntr


%  Outputs:
%  -cfg: configuration structure corrected and complemented
%  -method: 1,2,or 3 depending on whether cfg.method =  'trial', 
%           'trialTime', or 'ensemble' respectively 
%  -measInds: the indexes of the measures to be calculated
%  -Nmeasures: how many measures are to be calculated
%  -thisfunCommands: cell of the command strings for the calculation of
%                    the selected measures
%  -MeasNames: the cell with the names of all (sub)measures in the results'
% %             structure as constructed here:

% % -NmeasPmeas: the vector of the numbers of submeasures per measure 
% %              in the results' structure 


%Data validation, some unpacking and calculation of parameters an commands
funcName='DPstateCouplingPrepare';


if nargin<4
    z=[];
    cfg.Partial=0;
    Nvar = 2;
else
    cfg.Partial=1;
    Nvar = 3;
end


varName='cfg';
testCFG = {@(cfg)isstruct(cfg) };
param = {{}};
mode=['e'];
execfun={{}};
default=nan;
[RES, cfg] = DPvalidateData(cfg,testCFG,param,mode,execfun,default,varName,funcName);


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
    [RES, cfg.fs] = DPvalidateData(cfg.fs,testFS,param,mode,execfun,default,varName,funcName);
    cfg.Ts=1/cfg.fs; %sampling time in secs
else
    error('Structure field cfg.fs is missing.')
end

if isfield(cfg,'fc')
    varName='fc';
    testFC = {@(fc)isnumeric(fc),...
              @(fc)isreal(fc),...
              @(fc)all( fc>0 ),...
              @(fc)isvector(fc),...
              @(fc,Partial)( (Partial==0) && (length(fc)==2) ) || ( (Partial==1) && (length(fc)==3) ) };
    param={{},{},{},{},{cfg.Partial}};
    mode=['e','e','e','e','e'];
    execfun={{},{},{},{},{}};
    default=nan;
    [RES, cfg.fc] = DPvalidateData(cfg.fc,testFC,param,mode,execfun,default,varName,funcName);
    cfg.Tc=1./cfg.fc; %central period in secs
    cfg.TcMax = max(cfg.Tc); %and for the slower signal in case of CFC
    cfg.NtcMax = round(cfg.TcMax/cfg.Ts); %and in time points (samples)
else
    error('Structure field cfg.fc is missing.')
end

xSIZE = size(x);
xSIZElen = length(xSIZE);

varName='x';
testX = {@(x)isnumeric(x),...
         @(x)isreal(x),...
         @(x,xSIZElen) (xSIZElen==2)};

param={{},{},{xSIZElen}};
mode=['e','e','e'];
execfun={{},{},{}};
default=nan;
[RES, x] = DPvalidateData(x,testX,param,mode,execfun,default,varName,funcName);

cfg.N = xSIZE(1); %Number of time points
cfg.Ntr = xSIZE(2);%number of trials


ySIZE = size(y);
ySIZElen = length(ySIZE);
varName='y';
testY = {@(y)isnumeric(y),...
         @(y)isreal(y),...
         @(y,ySIZElen)(ySIZElen==2) ,...
         @(y,ySIZE,N,Ntr) all(ySIZE==[N Ntr]) };
param={{},{},{ySIZElen},{ySIZE,cfg.N,cfg.Ntr}};
mode=['e','e','e','e'];
execfun={{},{},{},{}};
default=nan;
[RES, y] = DPvalidateData(y,testY,param,mode,execfun,default,varName,funcName);

if (cfg.Partial==1)
    zSIZE = size(z);
    zSIZElen = length(zSIZE);
    varName='z';
    testZ = {@(z)isnumeric(z),...
             @(z)isreal(z),...
             @(z,zSIZElen)(zSIZElen==2) ,...
             @(z,zSIZE,N,Ntr) all(zSIZE==[N Ntr]) };
    param={{},{},{zSIZElen},{zSIZE,cfg.N,cfg.Ntr}};
    mode=['e','e','e','e'];
    execfun={{},{},{},{}};
    default=nan;
    [RES, z] = DPvalidateData(z,testZ,param,mode,execfun,default,varName,funcName);
end


if isfield(cfg,'time')
varName='time';
timSIZ = size(cfg.time);
testTIME = { @(time)isnumeric(time),...
             @(time)isreal(time),...
             @(time)isvector(time),... 
             @(time)all(diff(time)>0),...
             @(time)timSIZ(1)>=timSIZ(2),...
             @(time,N)length(time)<=N};
param={{},{},{},{},{},{cfg.N}};
mode=['e','e','e','e','w','w'];
execfun={{},{},{},{},@(time)time(:),@(time)time(1:cfg.N)};
default=([0:cfg.N-1]*cfg.Ts).';
[RES, cfg.time] = DPvalidateData(cfg.time,testTIME,param,mode,execfun,default,varName,funcName);
else
    %cprintf('Magenta','WARNING: Structure field cfg.time is missing. Setting default [0:N-1]/fs.\n')
    fprintf('WARNING: Structure field cfg.time is missing. Setting default [0:N-1]/fs.\n')
    cfg.time = ([0:cfg.N-1]*cfg.Ts).';
end
cfg.timeLen = cfg.time(end)-cfg.time(1);


if isfield(cfg,'D')
    varName='D';
    testD = {@(D)isnumeric(D),...
             @(D)isreal(D),...
             @(D)all( D>0 ),...
             @(D)isvector(D),...
             @(D,Nvar) (length(D)==Nvar),...
             @(D) all(round(D)==D)};
    param={{},{},{},{},{Nvar},{}};
    mode=['e','e','e','e','e','e'];
    execfun={{},{},{},{},{},{}};
    default=nan;
    [RES, cfg.D] = DPvalidateData(cfg.D,testD,param,mode,execfun,default,varName,funcName);
else
    error('Structure field cfg.D is missing.')
end


for iV = 1:Nvar;
    defaultTAU{iV}=ones(1,cfg.D(iV)-1)/cfg.fs;
end
if isfield(cfg,'tau')
    varName='tau';
    testTAU = {@(tau) iscell(tau),...
               @(tau,Nvar) (numel(tau)==Nvar),...
               @(tau) all(cellfun( @(tau_i)isnumeric(tau_i), tau )),...
               @(tau) all(cellfun( @(tau_i)isreal(tau_i), tau )),...
               @(tau,Ts) all(cellfun( @(tau_i)all(tau_i>=Ts), tau )),...
               @(tau) all(cellfun( @(tau_i)isvector(tau_i), tau )),...
               @(tau,D) all(cellfun( @(tau_i)all( numel(tau_i)==D-1 ), tau )),...
              };
    param={{},{Nvar},{},{},{cfg.Ts},{},{cfg.D}};
    mode=['e','e','e','e','e','e','e'];
    execfun={{},{},{},{},{},{},{}};
    [RES, cfg.tau] = DPvalidateData(cfg.tau,testTAU,param,mode,execfun,defaultTAU,varName,funcName);
else
    %cprintf('Magenta','WARNING: Structure field cfg.tau is missing. Setting all relative time delays equal to one sample point.\n')
    fprintf('WARNING: Structure field cfg.tau is missing. Setting all relative time delays equal to one sample point.\n')
    cfg.tau = defaultTAU;
end
for iV = 1:Nvar;
    cfg.tau{iV} = cumsum(cfg.tau{iV}); %cumulative embedding delays
    cfg.tauP{iV} = floor( cfg.tau{iV}*cfg.fs ); %...in sample points
    cfg.tau{iV} = cfg.tauP{iV}*cfg.Ts; %update tau in secs
    cfg.Nps(iV) = cfg.N-cfg.tauP{iV}(end);%...number of time points left for each variable after embedding
    cfg.maxTau(iV) = cfg.tau{iV}(end);
end

%Find number of embedded points after embedding of all time series
cfg.Np = min(cfg.Nps);
Np2=floor(cfg.Np/2); %a useful constant

%Calculate the point time vector
cfg.timeP = cfg.time(1:cfg.Np);

%Calculate the length of the point time:
cfg.timePlen = cfg.timeP(end)-cfg.timeP(1);

defaultWth= cfg.TcMax;
if isfield(cfg,'Wth')
    varName='Wth';
    testWth = {@(Wth)isnumeric(Wth),...
               @(Wth)isreal(Wth),...
               @(Wth)isscalar(Wth),...
               @(Wth)( Wth>=0 ) };
    param={{},{},{},{}};
    mode=['e','e','e','e'];
    execfun={{},{},{},{}};
    [RES, cfg.Wth] = DPvalidateData(cfg.Wth,testWth,param,mode,execfun,defaultWth,varName,funcName);
else
    cprint('Magenta','WARNING: Structure field cfg.Wth is missing. Setting default equal to the maximum main time scale (max(Tc)): %f secs.\n',defaultWth)
    cfg.Wth = defaultWth;
end
cfg.WthP = floor(cfg.Wth*cfg.fs);%Theiler window in sample points
cfg.Wth = cfg.WthP*cfg.Ts;%Theiler window in secs updated

if isfield(cfg,'normal')
    varName='normal';
    testNORMAL = { @(normal)ischar(normal),...
                   @(normal)isvector(normal),...
                   @(normal)strcmpi(normal,'none') ||strcmpi(normal,'zscore') || strcmpi(normal,'meanCntr') || strcmpi(normal,'linear') };
    param={{},{},{}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default=nan;
    [RES, cfg.normal] = DPvalidateData(cfg.normal,testNORMAL,param,mode,execfun,default,varName,funcName);
    
    if strcmpi(normal,'linear')
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
            [RES, cfg.normVal] = DPvalidateData(cfg.normVal,testNORMVAL,param,mode,execfun,default,varName,funcName);
        else
            cfg.normVal = [-1 1];
        end
    end
    
else
    cfg.normal = 'none';
end

if strcmpi(cfg.normal,'zscore')
    [x,cfg.xM,cfg.xSTD] = zscore(x);
    [y,cfg.yM,cfg.ySTD] = zscore(y);
    if (cfg.Partial==1)
        [z,cfg.zM,cfg.zSTD] = zscore(z);
    end
    
elseif strcmpi(cfg.normal,'meanCntr')
    cfg.xM = mean(x);
    x = x-repmat(cfg.xM,N,1);
    
    cfg.yM = mean(y);
    y = y-repmat(cfg.yM,N,1);
    
    if (cfg.Partial==1)
        z = z-repmat(cfg.zM,N,1);
    end
    
elseif strcmpi(cfg.normal,'linear')
    cfg.minX = min(x);
    cfg.maxX = max(x);
    range = cfg.maxX - cfg.minX;
    x = normVal(1) + repmat((diff(normalVal)./range),N,1) ./ ( x- repmat(cfg.minX,N,1) );
    
    cfg.minY = min(y);
    cfg.maxY = max(y);
    range = cfg.maxY - cfg.minY;
    y = normVal(1) + repmat((diff(normalVal)./range),N,1) ./ ( y- repmat(cfg.minY,N,1) );
    
    if (cfg.Partial==1)
        cfg.minZ = min(z);
        cfg.maxZ = max(z);
        range = cfg.maxZ - cfg.minZ;
        z = normVal(1) + repmat((diff(normalVal)./range),N,1) ./ ( z- repmat(cfg.minZ,N,1) );
    end
end


if isfield(cfg,'method')
    varName='method';
    testMETHOD = { @(method)ischar(method),...
                   @(method)isvector(method),...
                   @(method)strcmpi(method,'trial') || strcmpi(method,'trialTime') || strcmpi(method,'ensemble') };
    param={{},{},{}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default=nan;
    [RES, cfg.method] = DPvalidateData(cfg.method,testMETHOD,param,mode,execfun,default,varName,funcName);
else
    error('Structure field cfg.method is missing.')
end
method = 1*strcmpi(cfg.method,'trial') + 2*strcmpi(cfg.method,'trialTime') + 3*strcmpi(cfg.method,'ensemble');


if (method~=1)  
     
    if isfield(cfg,'timeCalc')
        varName='timeCalc';
        timCalcSIZ = size(cfg.timeCalc);
        testTIMECALC = { @(timeCalc)isnumeric(timeCalc),...
                         @(timeCalc)isreal(timeCalc),...
                         @(timeCalc)isvector(timeCalc),...
                         @(timeCalc)all(diff(timeCalc)>0),...
                         @(timeCalc,timCalcSIZ)timCalcSIZ(1)>=timCalcSIZ(2),...
                         @(timeCalc,N)length(timeCalc)<=N,...
                         @(timeCalc,time) all( ismember(timeCalc,time) )};
        param={{},{},{},{},{timCalcSIZ},{cfg.N},{cfg.time}};
        mode=['e','e','w','e','e','e','e'];
        execfun={{},{},@(time)time(:),{},{},{},{}};
        default=cfg.time;
        [RES, cfg.timeCalc] = DPvalidateData(cfg.timeCalc,testTIMECALC,param,mode,execfun,default,varName,funcName);
    else
        %cprintf('Magenta','WARNING: Structure field cfg.timeCalc is missing. Setting default cfg.timeCalc = cfg.timeP.\n')
        fprintf('WARNING: Structure field cfg.timeCalc is missing. Setting default cfg.timeCalc = cfg.timeP.\n')
        cfg.timeCalc = cfg.time;
    end
    if any(cfg.timeCalc>cfg.timeP(end))
        cfg.timeCalc = cfg.timeCalc(cfg.timeCalc<cfg.timeP(end));  %get rid of any excessive points after embedding
        cfg.timeCalc(end+1) = cfg.timeP(end);%...and set the last available point as a calculation point
    end
    cfg.Ncalc = length(cfg.timeCalc);
    
    %Initialize the time indices of calculation to be equal with all times:
    cfg.timeCalcIND = 1:cfg.Np;
    if ~isequal(cfg.timeCalc,cfg.timeP)
        if isfield(cfg,'upsample')
            varName='upsample';
            testUPSAMPLE = { @(upsample)ischar(upsample),...
                             @(upsample)isvector(upsample),...
                             @(upsample)strcmpi(upsample,'yes')|| strcmpi(upsample,'no') };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='no';
            [RES, cfg.upsample] = DPvalidateData(cfg.upsample,testUPSAMPLE,param,mode,execfun,default,varName,funcName);
        else
            cfg.upsample='no';
            %cprintf('Magenta','WARNING: Structure field cfg.upsample is missing. Setting default ''no''.\n')
            fprintf('WARNING: Structure field cfg.upsample is missing. Setting default ''no''.\n')
        end
        %Update cfg.timeCalcIND
        cfg.timeCalcIND = cfg.timeCalcIND( ismember(cfg.timeP,cfg.timeCalc) );
    else
        cfg.upsample='no';
    end
    
    if isfield(cfg,'winLen')
        varName='winLen';
        testWINLEN = { @(winLen)isnumeric(winLen),...
                       @(winLen)isreal(winLen),...
                       @(winLen)isscalar(winLen),...
                       @(winLen,timePLen)(winLen>=0)&&(winLen<timePLen)};
        param={{},{},{},{cfg.timePLen}};
        mode=['e','e','e','e'];
        execfun={{},{},{},{}};
        default=cfg.Ts;
        [RES, cfg.winLen] = DPvalidateData(cfg.winLen,testWINLEN,param,mode,execfun,default,varName,funcName);
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
    
    if isfield(cfg,'smoothWinfun')&&( isequal(cfg.timeCalc,cfg.timeP) || strcmpi(cfg.upsample,'yes') )
        varName='smoothWinfun';
        testSMOOTHWINFUN = { @(smoothWinfun)isa(smoothWinfun,'function_handle') };
        param={{}};
        mode=['e'];
        execfun={{}};
        default=@hanning;
        [RES, cfg.smoothWinfun] = DPvalidateData(cfg.smoothWinfun,testSMOOTHWINFUN,param,mode,execfun,default,varName,funcName);

        if isfield(cfg,'smoothWinlen')
            varName='smoothWinlen';
            testSMOOTHWINLEN = { @(smoothWinlen)isnumeric(smoothWinlen),...
                                 @(smoothWinlen)isreal(smoothWinlen),...
                                 @(smoothWinlen)isscalar(smoothWinlen),...
                                 @(smoothWinlen,Ts,timePLen)(smoothWinlen>Ts)&&(smoothWinlen<timePLen)};
            param={{},{},{},{cfg.Ts,cfg.timePLen}};
            mode=['e','e','e','e'];
            execfun={{},{},{},{}};
            default=min(cfg.Tc)/4;
            [RES, cfg.smoothWinlen] = DPvalidateData(cfg.smoothWinlen,testSMOOTHWINLEN,param,mode,execfun,default,varName,funcName);
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
                     @(measures)~isempty(intersect(measures,{'all','TE','MI','TI','NI','MP','CCD','PTE','PMI'})) };
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    default='all';
    [RES cfg.measures] = DPvalidateData(cfg.measures,testMEASURES,param,mode,execfun,default,varName,funcName);
    
else
    %cprintf('Magenta','WARNING: Structure field cfg.measures is missing. Replacing with ''all''.\n')
    fprintf('WARNING: Structure field cfg.measures is missing. Replacing with ''all''.\n')    
    cfg.measures='all';
end


%Cell of names of result's structure
MeasNames={'TE';...
'MI';...
'PTE';...
'PMI';...
'TI';...
'CCD';...
'SL';...
'NI';...
'MP';...
};
NmeasPmeas = [1 1 1 1 1 1 1 1 1]; %Number of 'measures per measures' categories

%Unpack the measures that will be calculated
if strcmpi(cfg.measures,'all')
    measures = ones(1,8);
else
    measures = [any(strcmpi(cfg.measures,'TE')),...
                any(strcmpi(cfg.measures,'MI')),...
				any(strcmpi(cfg.measures,'PTI')),...
				any(strcmpi(cfg.measures,'PMI'))
				any(strcmpi(cfg.measures,'TI')),...
				any(strcmpi(cfg.measures,'CCD')),...
				any(strcmpi(cfg.measures,'SL')),...
				any(strcmpi(cfg.measures,'NI')),...
				any(strcmpi(cfg.measures,'MP')),...
];
end

if any(measures([3,4])) && (cfg.Partial==0)
    warning('There is no third signal to calculate partial measures. They are ignored.')
    measures([3,4])=[0 0];
end

measInds = find(measures); %Indexes of measures to be calculated
Nmeasures = length(measInds); %number of measures to calculate
if (Nmeasures==0)
    error('Exiting because there are no measures to calculate!')
end


  
% %Measures' parameters structures:
% if measures(1)
%     if isfield(cfg,'IC')
%         varName='IC';
%         testIC = {@(IC)isstruct(IC) };
%         param = {{}};
%         mode=['e'];
%         execfun={{}};
%         default=nan;
%         [RES, cfg.IC] = DPvalidateData(cfg.IC,testIC,param,mode,execfun,default,varName,funcName);
%         
%         if isfield(cfg.IC,'Dphi0')
%             varName='IC.Dphi0';
%             testICDPHI0 = {@(Dphi0)isnumeric(Dphi0),...
%                            @(Dphi0)isreal(Dphi0),...
%                            @(Dphi0)( Dphi0>0 )&&( Dphi0<=pi/2 ),...
%                            @(Dphi0)isscalar(Dphi0) };
%             param={{},{},{},{}};
%             mode=['e','e','e','e'];
%             execfun={{},{},{},{}};
%             default=pi/4;
%             [RES, cfg.IC.Dphi0] = DPvalidateData(cfg.IC.Dphi0,testICDPHI0,param,mode,execfun,default,varName,funcName);
%         else
%             cprintf('Magenta','WARNING: Structure field cfg.IC.Dphi0 is missing. Replacing with default value pi/4.\n')
%             cfg.IC.Dphi0 = pi/4;
%         end
%     else
%         cprintf('Magenta','WARNING: Structure field cfg.IC is missing. Replacing with default value cfg.IC.Dphi0=pi/4.\n')
%         cfg.IC.Dphi0=pi/4;
%     end
% end





        
% %Calculate command strings
% if (method==1) %trial
%     
%     funCommands  = {'[PC.(''PCI'')(iTr), PC.(''NCI'')(iTr), PC.(''ACI'')(iTr), PC.(''ICI12'')(iTr), PC.(''ICI21'')(iTr)] = DPcalcIC(thisDphi,cfg.TcMax,cfg.IC.Dphi0);',...
%                     
%                      };
% 
% elseif (method==2) %trialTime 
%  
%      funCommands  = {'[PC.(''PCI'')(iT,iTr), PC.(''NCI'')(iT,iTr), PC.(''ACI'')(iT,iTr), PC.(''ICI12'')(iT,iTr), PC.(''ICI21'')(iT,iTr)] = DPcalcIC(thisDphi,cfg.TcMax,cfg.IC.Dphi0);',...
%                      
%                     };
% 
%                 
% else %ensemble 
%     
%     funCommands  = {'[PC.(''PCI'')(iT), PC.(''NCI'')(iT), PC.(''ACI'')(iT), PC.(''ICI12'')(iT), PC.(''ICI21'')(iT)] = DPcalcICensemble(thisDphi,cfg.NtcMax,cfg.IC.Dphi0,thisNwin,cfg.Ntr);',...
%                      
%                     };
% end        

% %Select the function commands of the measures to be calculated
% thisfunCommands = funCommands(measInds); 
% Ncomnds = length(thisfunCommands);


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
                               @(statMethod)strcmpi(statMethod,{'trialshuffling'})||strcmpi(statMethod,{'timeshuffling'})||strcmpi(statMethod,{'phaseshuffling'})||strcmpi(statMethod,{'phaserandom'}) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='trialshuffling';
            [RESULT, cfg.stats.method] = DPvalidateData(cfg.stats.method,testSTATMETHOD,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.stats.method is missing. Setting default value ''trialshuffling''.\n')
            fprintf('WARNING: Structure field cfg.stats.method is missing. Setting default value ''trialshuffling''.\n')
            cfg.stats.method = 'trialshuffling';
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
%         if strcmpi(cfg.stats.corrMultComp,'BONF')
%             switch method
%                 case 1
%                     defaultNperm = cfg.Ntr*defaultNperm;
%                 case 2
%                     defaultNperm =  cfg.Ntr*cfg.Ncalc*defaultNperm;
%                 case 3
%                     defaultNperm =  cfg.Ncalc*defaultNperm;
%             end
%         end
        
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


if ~isempty(cfg.stats)
    switch cfg.stats.method
        case 'timeshuffling'
            cfg.stats.surrfun = @(phi,N,D)DPsurrShufflTime(phi,N,D);
        case 'phaseshuffling'
            cfg.stats.surrfun = @(phi,N,D)DPsurrShufflPhase(phi,N,D);
        case 'phaserandom'
            cfg.stats.surrfun = @(phi,N,D)DPsurrRandPhase(phi,N,D); 
        otherwise 
            cfg.stats.surrfun = @(phi,N,D)DPsurrTrialShuf(phi,N,D);    
    end
end
