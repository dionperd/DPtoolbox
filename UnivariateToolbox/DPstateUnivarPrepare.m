function [cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, X] = DPstateUnivarPrepare(cfg,x)

%This function prepares data and configuration struture for DPstateUnivar
%

%Inputs: 
%--------------------------------------------------------------------------
%-cfg: configuration structure with fields
%   -fs: sampling frequency in herz, positive real scalar
%   -fc: central frequency of the signals in herz, 2 (or 3) element vector of
%        positive real values
%   -domain: 'real' or 'state', string, if domain=='real'
%   -D: number of time delay dimensions of x,
%   -tauX: a vector of D-1 embedding time delays of x, in samples
%   -time: time vector in secs, real valued vector, default=[0:N-1]/fs
%   -Wth: the Theiler window to be used for the nearest neighbors
%         statistics in samples
%--------------------------------------------------------------------------
%Optionally if domain=='real' then a time delay embedding is needed:
%        Then, Np is the number of points in the pointset
%--------------------------------------------------------------------------
%   -normal: 'zscore' or 'meanCenter', 'linear', string, default='none' 
%           if normal == 'linear', there should also be a vector of 2 real
%           numbers, normVal for the minimum and maximum value, 
%           default, normVal = [-1 1]
%
%           but only if domain='real' 
%--------------------------------------------------------------------------
%   -method: one of 
%            1. 'trial', for estimation per trial
%            2. 'trialTime', for estimation per trial with a sliding time 
%               window
%            3. 'ensemble', for estimation across trials, either pointwise
%            or with a time window

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
%              1.'CD' for correlation dimension 
%
%   -CD: optional structure of inputs related to CD
%       -k: number of nearest neighbors
%--------------------------------------------------------------------------
%Optionally if statistics are to be made:
%   -stats: a structure with fields:
%       -method:a string with the surrogate test method that is to be applied:
%                 -"phaseshuffling": random phase shuffling in Fourier
%                 space
%                default: "phaseshuffling"
%       -surrfun: a handle to a function of the form
%                 xSurr = surrfun(x,N,D), where
%                   -x: the data series columnwise matrix, 
%                   -N: the number of data points per dimension of the set,
%                    here it is N=cfg.Np
%                   -D: the number of dimensions of the set, here it is D=1
%                that creates a surrogate data set given the method above
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
% -x: the first signal, time series of N time points, 
%     either real valued signal: 
%     a matrix N x 1 x Ntr
%     or state space embedded points:
%     a matrix N x D x Ntr 
%     where D is the number of embedding dimensions, and Ntr the number of
%     trials



%  Outputs:
%  -cfg: configuration structure corrected and complemented
%  -method: 1,2,or 3 depending on whether cfg.method =  'trial', 
%           'trialTime', or 'ensemble' respectively 
%  -measInds: the indexes of the measures to be calculated
%  -Nmeasures: how many measures are to be calculated
%  -thisfunCommands: cell of the command strings for the calculation of
%                    the selected measures
%  -X: embedded pointset of signal x, matrix (Np,D,Ntr)




%Data validation, some unpacking and calculation of parameters and commands
funcName='DPstateUnivarPrepare';


varName='cfg';
testCFG = {@(cfg)isstruct(cfg) };
param = {{}};
mode=['e'];
execfun={{}};
default=nan;
[~, cfg] = DPvalidateData(cfg,testCFG,param,mode,execfun,default,varName,funcName);

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
    [~, cfg.fs] = DPvalidateData(cfg.fs,testFS,param,mode,execfun,default,varName,funcName);
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
              @(fc,MultiVar)( (MultiVar==0) && (length(fc)==2) ) || ( (MultiVar==1) && (length(fc)==3) ) };
    param={{},{},{},{},{cfg.MultiVar}};
    mode=['e','e','e','e','e'];
    execfun={{},{},{},{},{}};
    default=nan;
    [~, cfg.fc] = DPvalidateData(cfg.fc,testFC,param,mode,execfun,default,varName,funcName);
    cfg.Tc=1./cfg.fc; %central period in secs
    cfg.TcMax = max(cfg.Tc); %and for the slower signal in case of CFC
    cfg.NtcMax = round(cfg.TcMax/cfg.Ts); %and in time points (samples)
else
    error('Structure field cfg.fc is missing.')
end

if isfield(cfg,'domain')
    varName='domain';
    testDOMAIN = { @(domain)ischar(domain),...
                   @(domain)isvector(domain),...
                   @(domain)strcmpi(domain,'state')|| strcmpi(domain,'real') };
    param={{},{},{}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default=nan;
    [~, cfg.domain] = DPvalidateData(cfg.domain,testDOMAIN,param,mode,execfun,default,varName,funcName);
else
    error('Structure field cfg.domain is missing.')
end


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
[~, x] = DPvalidateData(x,testX,param,mode,execfun,default,varName,funcName);
N = xSIZE(1); %Number of time points
cfg.N =N;
Np=N; %initialize Np as equal to N
Np2 = floor(Np/2); %a usefull constant
cfg.Np=Np;
%Calculate the cut time vector
cfg.timeP = cfg.time;
%Calculate the duration of the cut time:
cfg.timePlen = cfg.timeP(end)-cfg.timeP(1);
if xSIZElen<3
    Ntr=xSIZE(2);
    D = 1;
else
    Ntr = xSIZE(3); %Number of trials
    D=xSIZE(2);
end
cfg.Ntr =Ntr;
cfg.D=D;

if strcmpi(cfg.domain,'real')
    if isfield(cfg,'normal')
        varName='normal';
        testNORMAL = { @(normal)ischar(normal),...
            @(normal)isvector(normal),...
            @(normal)strcmpi(normal,'none') ||strcmpi(normal,'zscore') || strcmpi(normal,'meanCntr') || strcmpi(normal,'linear') };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default=nan;
        [~, cfg.normal] = DPvalidateData(cfg.normal,testNORMAL,param,mode,execfun,default,varName,funcName);
        
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
                [~, cfg.normVal] = DPvalidateData(cfg.normVal,testNORMVAL,param,mode,execfun,default,varName,funcName);
            else
                cfg.normVal = [-1 1];
            end
        end
        
    else
        cfg.normal = 'none';
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
        x = normVal(1) + repmat((diff(normalVal)./range),N,1) ./ ( x- repmat(cfg.minX,N,1) );
        
    end
end


if isfield(cfg,'measures')
    varName='measures';
    testMEASURES = { @(measures)iscellstr(measures)||ischar(measures),...
                     @(measures)~isempty(intersect(measures,{'all','CD'})) };
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    default='all';
    [~, cfg.measures] = DPvalidateData(cfg.measures,testMEASURES,param,mode,execfun,default,varName,funcName);
    
else
    cprintf('Magenta','WARNING: Structure field cfg.measures is missing. Replacing with ''all''.\n')
    cfg.measures='all';
end


%Unpack the measures that will be calculated
if strcmpi(cfg.measures,'all')
    measures = ones(1,9);
else
    measures = [any(strcmpi(cfg.measures,'CD')),...
                 ];
end

measInds = find(measures); %Indexes of measures to be calculated
Nmeasures = length(measInds); %number of measures to calculate
if (Nmeasures==0)
    error('Exiting because there are no measures to calculate!')
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
param={{},{},{},{},{},{N}};
mode=['e','e','e','e','w','w'];
execfun={{},{},{},{},@(time)time(:),@(time)time(1:N)};
default=([0:N-1]*cfg.Ts).';
[~, cfg.time] = DPvalidateData(cfg.time,testTIME,param,mode,execfun,default,varName,funcName);
else
    cprintf('Magenta','WARNING: Structure field cfg.time is missing. Setting default [0:N-1]/fs.\n')
    cfg.time = ([0:N-1]*cfg.Ts).';
end
cfg.timeLen = cfg.time(end)-cfg.time(1);


% 
% if strcmpi(cfg.domain,'real')
%     
%     cfg.Np=Np;
%     %Calculate the point time vector
%     cfg.timeP = cfg.time(1:Np); 
%     %Calculate the duration of the point time:
%     cfg.timePlen = cfg.timeP(end)-cfg.timeP(1);
% end


if isfield(cfg,'method')
    varName='method';
    testMETHOD = { @(method)ischar(method),...
                   @(method)isvector(method),...
                   @(method)strcmpi(method,'trial') || strcmpi(method,'trialTime') || strcmpi(method,'ensemble') };
    param={{},{},{}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default=nan;
    [~, cfg.method] = DPvalidateData(cfg.method,testMETHOD,param,mode,execfun,default,varName,funcName);
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
                         @(timeCalc,Ncut)length(timeCalc)<=Ncut,...
                         @(timeCalc,timeCut) all( ismember(timeCalc,timeCut) )};
        param={{},{},{},{},{timCalcSIZ},{cfg.Ncut},{cfg.timeCut}};
        mode=['e','e','w','e','e','e','e'];
        execfun={{},{},@(time)time(:),{},{},{},{}};
        default=cfg.timeCut;
        [~, cfg.timeCalc] = DPvalidateData(cfg.timeCalc,testTIMECALC,param,mode,execfun,default,varName,funcName);
    else
        cprintf('Magenta','WARNING: Structure field cfg.timeCalc is missing. Setting default cfg.timeCalc = cfg.timeCut.\n')
        cfg.timeCalc = cfg.timeCut;
    end
    cfg.Ncalc = length(cfg.timeCalc);
    
    %Initialize the time indices of calculation to be equal with all times:
    cfg.timeCalcIND = 1:cfg.Ncut;
    if ~isequal(cfg.timeCalc,cfg.timeCut)
        if isfield(cfg,'upsample')
            varName='upsample';
            testUPSAMPLE = { @(upsample)ischar(upsample),...
                             @(upsample)isvector(upsample),...
                             @(upsample)strcmpi(upsample,'yes')|| strcmpi(upsample,'no') };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='no';
            [~, cfg.upsample] = DPvalidateData(cfg.upsample,testUPSAMPLE,param,mode,execfun,default,varName,funcName);
        else
            cfg.upsample='no';
            cprintf('Magenta','WARNING: Structure field cfg.upsample is missing. Setting default ''no''.\n')
        end
        %Update cfg.timeCalcIND
        cfg.timeCalcIND = cfg.timeCalcIND( ismember(cfg.timeCut,cfg.timeCalc) );
    else
        cfg.upsample='no';
    end
    
    if isfield(cfg,'winLen')
        varName='winLen';
        testWINLEN = { @(winLen)isnumeric(winLen),...
                       @(winLen)isreal(winLen),...
                       @(winLen)isscalar(winLen),...
                       @(winLen,timeCutLen)(winLen>=0)&&(winLen<timeCutLen)};
        param={{},{},{},{cfg.timeCutLen}};
        mode=['e','e','e','e'];
        execfun={{},{},{},{}};
        default=cfg.Ts;
        [~, cfg.winLen] = DPvalidateData(cfg.winLen,testWINLEN,param,mode,execfun,default,varName,funcName);
    else
        cprintf('Magenta','WARNING: Structure field cfg.winLen is missing. Setting cfg.winLen=Ts=%f.\n',cfg.Ts)
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
    
    if isfield(cfg,'smoothWinfun')&&( isequal(cfg.timeCalc,cfg.timeCut) || strcmpi(cfg.upsample,'yes') )
        varName='smoothWinfun';
        testSMOOTHWINFUN = { @(smoothWinfun)isa(smoothWinfun,'function_handle') };
        param={{}};
        mode=['e'];
        execfun={{}};
        default=@hanning;
        [~, cfg.smoothWinfun] = DPvalidateData(cfg.smoothWinfun,testSMOOTHWINFUN,param,mode,execfun,default,varName,funcName);

        if isfield(cfg,'smoothWinlen')
            varName='smoothWinlen';
            testSMOOTHWINLEN = { @(smoothWinlen)isnumeric(smoothWinlen),...
                                 @(smoothWinlen)isreal(smoothWinlen),...
                                 @(smoothWinlen)isscalar(smoothWinlen),...
                                 @(smoothWinlen,Ts,timeCutLen)(smoothWinlen>Ts)&&(smoothWinlen<timeCutLen)};
            param={{},{},{},{cfg.Ts,cfg.timeCutLen}};
            mode=['e','e','e','e'];
            execfun={{},{},{},{}};
            default=min(cfg.Tc)/4;
            [~, cfg.smoothWinlen] = DPvalidateData(cfg.smoothWinlen,testSMOOTHWINLEN,param,mode,execfun,default,varName,funcName);
        else
            cprintf('Magenta','WARNING: Structure field cfg.smoothWinlen is missing. Setting cfg.smoothWinlen=min(cfg.Tc)/4=%f.\n',min(cfg.Tc)/4)
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


    

%Measures' parameters structures:
if measures(1)
    if isfield(cfg,'CD')
        varName='CD';
        testCD = {@(CD)isstruct(CD) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [~, cfg.CD] = DPvalidateData(cfg.CD,testCD,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.CD,'k')
            varName='IC.k';
            testK = {@(k)isnumeric(k),...
                           @(k)isreal(k),...
                           @(k)( k>0 ),...
                           @(k)isscalar(k) };
            param={{},{},{},{}};
            mode=['e','e','e','e'];
            execfun={{},{},{},{}};
            default=pi/4;
            [~, cfg.IC.Dphi0] = DPvalidateData(cfg.IC.Dphi0,testICDPHI0,param,mode,execfun,default,varName,funcName);
        else
            cprintf('Magenta','WARNING: Structure field cfg.IC.Dphi0 is missing. Replacing with default value pi/4.\n')
            cfg.IC.Dphi0 = pi/4;
        end
    else
        cprintf('Magenta','WARNING: Structure field cfg.IC is missing. Replacing with default value cfg.IC.Dphi0=pi/4.\n')
        cfg.IC.Dphi0=pi/4;
    end
end



%Surrogate statistics structure:
if isfield(cfg,'stats')
    
    varName='stats';
    testSTATS = { @(stats)isstruct(stats)||isempty(stats) };
    param = {{}};
    mode=['e'];
    execfun={{}};
    default=nan;
    [~, cfg.stats] = DPvalidateData(cfg.stats,testSTATS,param,mode,execfun,default,varName,funcName);
    
    
    if ~isempty(cfg.stats)
        
        if isfield(cfg.stats,'method')
            varName='stats.method';
            testSTATMETHOD = { @(statMethod)ischar(statMethod),...
                @(statMethod)isvector(statMethod),...
                @(statMethod)strcmpi(statMethod,{'phaseshuffling'}) };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='phaseshuffling';
            [~, cfg.stats.method] = DPvalidateData(cfg.stats.method,testSTATMETHOD,param,mode,execfun,default,varName,funcName);
        else
            cprintf('Magenta','WARNING: Structure field cfg.stats.method is missing. Setting default value ''trialshuffling''.\n')
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
            default=0.01;
            [~, cfg.stats.alpha] = DPvalidateData(cfg.stats.alpha,testSTATSALPHA,param,mode,execfun,default,varName,funcName);
        else
            cprintf('Magenta','WARNING: Structure field cfg.stats.alpha is missing. Setting default value 0.01.\n')
            cfg.stats.alpha = 0.01;
        end
        
        if isfield(cfg.stats,'tail')
            varName='stats.tail';
            testSTATSTAIL = {@(tail) (tail==1) || (tail==2)};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=1;
            [~, cfg.stats.tail] = DPvalidateData(cfg.stats.tail,testSTATSTAIL,param,mode,execfun,default,varName,funcName);
        else
            cprintf('Magenta','WARNING: Structure field cfg.stats.tail is missing. Setting default value 1.\n')
            cfg.stats.tail = 1;
        end
        
        defaultNperm = cfg.stats.tail*ceil(1/cfg.stats.alpha)+1;
        if isfield(cfg.stats,'Nperm')
            varName='stats.Nperm';
            testSTATSNPERM = {@(Nperm)isnumeric(Nperm),...
                @(Nperm)isreal(Nperm),...
                @(Nperm)isscalar(Nperm) ,...
                @(Nperm)ceil(Nperm)==Nperm,...
                @(Nperm,defaultNperm)( Nperm>defaultNperm ) };
            
            param={{},{},{},{},{defaultNperm}};
            mode=['e','e','e','w','e'];
            execfun={{},{},{},@(Nperm)ceil(Nperm),{}};
            default=defaultNperm;
            [~, cfg.stats.Nperm] = DPvalidateData(cfg.stats.Nperm,testSTATSNPERM,param,mode,execfun,default,varName,funcName);
        else
            cprintf('Magenta','WARNING: Structure field cfg.stats.Nperm is missing. Setting default value tail*ceil(1/alpha)+1=%f.\n',defaultNperm)
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
            [~, cfg.stats.pointStatMethod] = DPvalidateData(cfg.stats.pointStatMethod,testPOINTSTATMETHOD,param,mode,execfun,default,varName,funcName);
        else
            cprintf('Magenta','WARNING: Structure field cfg.stats.pointStatMethod is missing. Setting default value ''''.\n')
            cfg.stats.pointStatMethod = '';
        end
        
        
        if isfield(cfg.stats,'multiStatfun')
            
            varName='stats.multiStatfun';
            testMULTISTATSFUN = {@(multiStatfun) isstruct(multiStatfun)||isempty(multiStatfun)};
            param={{}};
            mode=['e'];
            execfun={{}};
            default=[];
            [~, cfg.stats.multiStatfun] = DPvalidateData(cfg.stats.multiStatfun,testMULTISTATSFUN,param,mode,execfun,default,varName,funcName);
            
            if ~isempty(cfg.stats.multiStatfun)
                
                if isfield(cfg.stats.multiStatfun,'fun')
                    varName='stats.multiStatfun.fun';
                    testMULTISTATSFUNCTION = {@(fun) isa(fun,'function_handle')};
                    param={{}};
                    mode=['e'];
                    execfun={{}};
                    default=nan;
                    [~, cfg.stats.multiStatfun.fun] = DPvalidateData(cfg.stats.multiStatfun.fun,testMULTISTATSFUNCTION,param,mode,execfun,default,varName,funcName);
                    
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
        
        
        if isfield(cfg.stats,'corrMultComp')
            varName='stats.corrMultComp';
            testSTATMULTCOMP = {@(corrMultComp)ischar(corrMultComp),...
                @(corrMultComp)isvector(corrMultComp),...
                @(corrMultComp)any( [strcmpi(corrMultComp,{'BONF'}), strcmpi(corrMultComp,{'FDR'}), strcmpi(corrMultComp,{''})]  )  };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default={''};
            [~, cfg.stats.corrMultComp] = DPvalidateData(cfg.stats.corrMultComp,testSTATMULTCOMP,param,mode,execfun,default,varName,funcName);
        else
            cfg.stats.corrMultComp='';
        end
        
        
    end
    
else
    cfg.stats=[];
end

        
%Calculate command strings
if (method==1) %trial
    
    funCommands  = {'PC.(''CD'')(iTr) = DPcalcCD(thisX,cfg.k);', ...
                    };

elseif (method==2) %trialTime 
 
     funCommands  = {'PC.(''CD'')(iT,iTr) = DPcalcCD(thisX,cfg.k);',...
                     };

                
else %ensemble 
    
    funCommands  = {'PC.(''CD'')(iT) = DPcalcCDens(thisX,cfg.k);', ...
                    };
end   

%Select the function commands of the measures to be calculated
thisfunCommands = funCommands(measInds); 
Ncomnds = length(thisfunCommands);

if ~isempty(cfg.stats)
    %Create a matrix of handles to function of surrogate time series' creation:
    switch cfg.stats.method
        case 'phaseshuffling'
            cfg.stats.surrfun = @(Data,N,D)DPphaeShufSurr(Data,N,D);
    end
end
