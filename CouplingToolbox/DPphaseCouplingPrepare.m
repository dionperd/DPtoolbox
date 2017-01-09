function [phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas] = DPphaseCouplingPrepare(x,cfg)

%This function prepares data and configuration struture for DPphaseCoupling
%

%Inputs:

% -x: the time series of N time points, 
%     a matrix either Nx2 (for a single trial) or Nx2xNtr of Ntr trials data,
%     either real valued signal, or phase signal 

%-cfg: configuration structure with fields
%   -fs: sampling frequency in herz, positive real scalar
%   -fc: central frequency of the signals in herz, 2 element vector of
%        positive real values
%   -time: time vector in secs, real valued vector, default=[0:N-1]/fs
%   -domain: 'real' or 'phase', string, if domain=='real', they have to be
%            mean centered
%--------------------------------------------------------------------------
%Optionally if domain=='real' and a hilbert transform is needed:
%   -cutTails: parts of the time series at the begining and end to be cut 
%                   out of the calculation in secs, vector of 2
%                   0=<positive real<N/fs/2 values, default: [TcMax TcMax]
%                   where TcMax=1/min(fc)

%        Then, Ncut = N-round(cfg.cutBegin*cfg.fs)-round(cfg.cutEnd*cfg.fs)
%--------------------------------------------------------------------------
%   -nm of n:m coupling, matrix of dimension K x 2 of 
%              for K different calculations of phase coupling
%              positive integer numbers, default: lcm(cfg.fc(1),cfg.fc(2))/cfg.fc
%   -method: one of 
%            1. 'trial', for estimation per trial
%            2. 'trialTime', for estimation per trial with a sliding time 
%               window
%            3. 'ensemble', for estimation across trials, either pointwise
%            or with a time window
%   -phaseTrnsfrm: 'yes' or 'no',string, default='no' 
%                but only if (Ncut or winlen)>TcMax (look below) 
%   -Nbins: number of bins, positive real odd integer, 
%           default: calculated from data  
%           in the interval [5 33] as 2*floor(sqrt(N)/2)+1
%--------------------------------------------------------------------------
%Optionally for 'trialTime' or 'ensemble'
%   -timeCalc: a subset of timeCut (in secs), with the time points where the calculation
%       will be performed, default=timeCut
%   -winLen: a positive number for the time length of the time window of
%   calculation, default = 1/fs 
% leading to a window length of Nwin points
%
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
%--------------------------------------------------------------------------
%   -measures: the coupling measures that should be calculated
%              among the available ones:
%              1.'IC' for integrative coupling index
%                 and the rest of Viktor Mueller measures (PCI, NCI, ACI, ICI12, ICI21)
%                 with clearing period:
%                   -min([T,Ncut]) for method: 'trial'
%                   -min([T,Nwin]) for method: 'trialTime' and 'ensemble'
%                 a warning is printed in case that clearing
%                 period<TcMax=1/min(fc)
%              2.'PLV' for phase locking value/index
%              3.'PPC' for pairwise phase consistency
%              4.'PLI' for phase lag index
%              5.'SE' for index based on Shannon entropy,
%              6.'CP' for index based on conditional probability,
%              7.'MI' for mutual information
%              8.'PR' for index based on phase reconstruction (only with
%                     as long as (Ncut or winlen)>TcMax=1/min(fc)
%              9.'CCR' for circular correlation coefficient

%                 cell of strings, default {'IC','PLV','PPC','PLI','SE','CP','MI''PR','CCR'} or
%                 string 'all' for all of them, default: 'all'

%   -IC: optional structure of inputs related to IC
%       -Dphi0: phase sync threshold, 0<positive real<pi/2, default=pi/4
%   -PR: optional structure of inputs related to PR
%       -method: string either 'fourier' or 'iter', default: 'fourier'
%       -Ngrid: number of grid points, positive real integer,
%               default: calculated from data
%       -or: order of fourier expansion, positive real, default=10
%--------------------------------------------------------------------------
%Optionally if statistics are to be made:
%   -stats: a structure with fields:
%       -method:a string with the surrogate test method that is to be applied:
%                 -"trialshuffling": random trial shuffling
%                 -"timeshuffling": random time shuffling 
%                 -"phaseshuffling": random phase shuffling in Fourier
%                 space
%                 -"phaserandom": randomized phase in Fourier space
%                default: "trialshuffling"
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
%                                 multiStat = multiStatfun(PC,PCsurr,pointStat,params) where
%                                 -PC: the result of the original data
%                                 -PCsurr: the results of the surrogate data
%                                 -pointStat: a structure of the form of PC
%                                             with values of the point statistic
%                                 -multiStat: a structure of the form of PC
%                                             with values a multivariate statistic
%                     -params: parameters for that function



%  Outputs:

%  -phi: phases of the signals wrapped into the interval [0 2pi)
%  -cfg: configuration structure corrected and complemented
%  -method: 1,2,3 or 4 depending on whether cfg.method =  'trial', 
%           'trialWin','time', or 'timeWin' respectively 
%  -measInds: the indexes of the measures to be calculated
%  -Nmeasures: how many measures are to be calculated
%  -thisfunCommands: cell of the command strings for the calculation of
%                    the selected measures
%  -MeasNames: the cell with the names of all (sub)measures in the results'
% %             structure as constructed here:

% % -NmeasPmeas: the vector of the numbers of submeasures per measure 
% %              in the results' structure 



%Define a constant:
TwoPi = 2*pi;


%Data validation, some unpacking and calculation of parameters an commands
funcName='DPphaseCouplingPrepare';

xSIZE = size(x);
xSIZElen = length(xSIZE);
varName='x';
testX = {@(x)isnumeric(x),...
         @(x)isreal(x),...
         @(x,xSIZElen)xSIZElen<4,...
         @(x,xSIZE)xSIZE(2)==2 };

param={{},{},{xSIZElen},{xSIZE}};
mode=['e','e','e','e'];
execfun={{},{},{},{}};
default=nan;
[RESULT, x] = DPvalidateData(x,testX,param,mode,execfun,default,varName,funcName);
N = xSIZE(1); %Number of time points
cfg.N =N;
Ncut=N; %initialize Ncut as equal to N
Ncut2 = floor(Ncut/2); %a usefull constant
cfg.Ncut=Ncut;
%Calculate the cut time vector
cfg.timeCut = cfg.time;
%Calculate the duration of the cut time:
cfg.timeCutLen = cfg.timeCut(end)-cfg.timeCut(1);
if xSIZElen<3
    Ntr=1;
else
    Ntr = xSIZE(3); %Number of trials
end
cfg.Ntr =Ntr;
% Xmax = max(x(:));
% Xmin = min(x(:));


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

if isfield(cfg,'fc')
    varName='fc';
    testFC = {@(fc)isnumeric(fc),...
              @(fc)isreal(fc),...
              @(fc)all( fc>0 ),...
              @(fc)isvector(fc),...
              @(fc)length(fc)==2 };
    param={{},{},{},{},{}};
    mode=['e','e','e','e','e'];
    execfun={{},{},{},{},{}};
    default=nan;
    [RESULT, cfg.fc] = DPvalidateData(cfg.fc,testFC,param,mode,execfun,default,varName,funcName);
    cfg.Tc=1./cfg.fc; %central period in secs
    cfg.TcMax = max(cfg.Tc); %and for the slower signal in case of CFC
    cfg.NtcMax = round(cfg.TcMax/cfg.Ts); %and in time points (samples)
else
    error('Structure field cfg.fc is missing.')
end


if isfield(cfg,'time')
varName='time';
testTIME = { @(time)isnumeric(time),...
             @(time)isreal(time),...
             @(time)isvector(time),... 
             @(time)all(diff(time)>0),...
             @(time,N)length(time)<=N};
param={{},{},{},{},{N}};
mode=['e','e','w','e','w'];
execfun={{},{},@(time)time(:),{},@(time)time(1:N)};
default=([0:N-1]*cfg.Ts).';
[RESULT, cfg.time] = DPvalidateData(cfg.time,testTIME,param,mode,execfun,default,varName,funcName);
else
    %cprintf('Magenta','WARNING: Structure field cfg.time is missing. Setting default [0:N-1]/fs.\n')
    fprintf('WARNING: Structure field cfg.time is missing. Setting default [0:N-1]/fs.\n')
    cfg.time = ([0:N-1]*cfg.Ts).';
end
cfg.timeLen = cfg.time(end)-cfg.time(1);

if isfield(cfg,'domain')
    varName='domain';
    testDOMAIN = { @(domain)ischar(domain),...
                   @(domain)isvector(domain),...
                   @(domain)strcmpi(domain,'phase')|| strcmpi(domain,'real') };
    param={{},{},{}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default=nan;
    [RESULT, cfg.domain] = DPvalidateData(cfg.domain,testDOMAIN,param,mode,execfun,default,varName,funcName);
else
    error('Structure field cfg.domain is missing.')
end


if strcmpi(cfg.domain,'real')
    
    if isfield(cfg,'cutTails')
        varName='cutTails';
        testCUTTAILS = {@(cutTails)isnumeric(cutTails),...
                        @(cutTails)isreal(cutTails),...
                        @(cutTails,timeLen)all( cutTails>=0 & cutTails<timeLen/2 ),...
                        @(cutTails)isvector(cutTails),...
                        @(cutTails)length(cutTails)==2 };
        param={{},{},{cfg.timeLen},{},{}};
        mode=['e','e','e','e','e'];
        execfun={{},{},{},{},{}};
        default=[cfg.TcMax cfg.TcMax];
        [RESULT, cfg.cutTails] = DPvalidateData(cfg.cutTails,testCUTTAILS,param,mode,execfun,default,varName,funcName);
    else
        %cprintf('Magenta','WARNING: Structure field cfg.cutTails is missing. Setting default [1/min(fc) 1/min(fc)].\n')
        fprintf('WARNING: Structure field cfg.cutTails is missing. Setting default [1/min(fc) 1/min(fc)].\n')
        cfg.cutTails=[cfg.TcMax cfg.TcMax];
    end

    %Update Ncut after cutting
    cfg.cutTailsPoints = round(cfg.cutTails*cfg.fs);
    Ncut = N-sum(cfg.cutTailsPoints);
    Ncut2=floor(Ncut/2);
    cfg.Ncut=Ncut;
    %Calculate the cut time vector
    cfg.timeCut = cfg.time(cfg.cutTailsPoints(1)+1:N-cfg.cutTailsPoints(2)); 
    %Calculate the duration of the cut time:
    cfg.timeCutLen = cfg.timeCut(end)-cfg.timeCut(1);
end


if isfield(cfg,'nm')
varName='nm';
nmSIZE=size(cfg.nm);
testNM = { @(nm)isnumeric(nm),...
           @(nm)isreal(nm),...
           @(nm)all(nm(:)>0),...
           @(nm)all(round(nm(:))==nm(:)),...
           @(nm,nmSIZE)length(nmSIZE)==2,...
           @(nm,nmSIZE)nmSIZE(2)==2};
param={{},{},{},{},{nmSIZE},{nmSIZE}};
mode=['e','e','e','e','e','e'];
execfun={{},{},{},{},{},{}};
fcint = floor(10^6*cfg.fc); %integer fc to a precision of 10^(-6)
default=fcint([2 1])/gcd( fcint(1), fcint(2) );
[RESULT, cfg.nm] = DPvalidateData(cfg.nm,testNM,param,mode,execfun,default,varName,funcName);

% Find the unique rows of cfg.nm
cfg.nm = unique(cfg.nm, 'rows');

else
    %cprintf('Magenta','WARNING: Structure field cfg.nm is missing. Setting default [1 1].\n')
    fprintf('WARNING: Structure field cfg.nm is missing. Setting default [1 1].\n')
    cfg.nm = [1 1];
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
    [RESULT, cfg.method] = DPvalidateData(cfg.method,testMETHOD,param,mode,execfun,default,varName,funcName);
else
    error('Structure field cfg.method is missing.')
end
method = 1*strcmpi(cfg.method,'trial') + 2*strcmpi(cfg.method,'trialTime') + 3*strcmpi(cfg.method,'ensemble');


if (method~=1)  
    
    
    if isfield(cfg,'timeCalc')
        varName='timeCalc';
        testTIMECALC = { @(timeCalc)isnumeric(timeCalc),...
                         @(timeCalc)isreal(timeCalc),...
                         @(timeCalc)isvector(timeCalc),...
                         @(timeCalc)all(diff(timeCalc)>0),...
                         @(timeCalc,Ncut)length(timeCalc)<=Ncut,...
                         @(timeCalc,timeCut) all( ismember(timeCalc,timeCut) )};
        param={{},{},{},{},{cfg.Ncut},{cfg.timeCut}};
        mode=['e','e','w','e','e','e'];
        execfun={{},{},@(time)time(:),{},{},{}};
        default=cfg.timeCut;
        [RESULT, cfg.timeCalc] = DPvalidateData(cfg.timeCalc,testTIMECALC,param,mode,execfun,default,varName,funcName);
    else
        %cprintf('Magenta','WARNING: Structure field cfg.timeCalc is missing. Setting default cfg.timeCalc = cfg.timeCut.\n')
        fprintf('WARNING: Structure field cfg.timeCalc is missing. Setting default cfg.timeCalc = cfg.timeCut.\n')
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
            [RESULT, cfg.upsample] = DPvalidateData(cfg.upsample,testUPSAMPLE,param,mode,execfun,default,varName,funcName);
        else
            cfg.upsample='no';
            %cprintf('Magenta','WARNING: Structure field cfg.upsample is missing. Setting default ''no''.\n')
            fprintf('WARNING: Structure field cfg.upsample is missing. Setting default ''no''.\n')
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
    %This is the time window to check for the protophase transformation and
    %the PR measures
    TPR = cfg.winLen;
    
    if isfield(cfg,'smoothWinfun')&&( isequal(cfg.timeCalc,cfg.timeCut) || strcmpi(cfg.upsample,'yes') )
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
                                 @(smoothWinlen,Ts,timeCutLen)(smoothWinlen>Ts)&&(smoothWinlen<timeCutLen)};
            param={{},{},{},{cfg.Ts,cfg.timeCutLen}};
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
else
    %This is the time window to check for the protophase transformation and
    %the PR measures
    TPR = cfg.timeCutLen;
end

%Make sure that the window is at least one period long for protophases to
%be transformed to true phases
if (TPR>=cfg.TcMax)
    if isfield(cfg,'phaseTrnsfrm')
        varName='phaseTrnsfrm';
        testPHASETRNSFRM = { @(phaseTrnsfrm)ischar(phaseTrnsfrm),...
                             @(phaseTrnsfrm)isvector(phaseTrnsfrm),...
                             @(phaseTrnsfrm)strcmpi(phaseTrnsfrm,'yes')|| strcmpi(phaseTrnsfrm,'no') };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='no';
        [RESULT, cfg.phaseTrnsfrm] = DPvalidateData(cfg.phaseTrnsfrm,testPHASETRNSFRM,param,mode,execfun,default,varName,funcName);
    else
        cfg.phaseTrnsfrm='no';
        %cprintf('Magenta','WARNING: Structure field cfg.phaseTrnsfrm is missing. Setting default ''no''.\n')
        fprintf('WARNING: Structure field cfg.phaseTrnsfrm is missing. Setting default ''no''.\n')
    end
else
    cfg.phaseTrnsfrm='no';
end
cfg.phaseTrnsfrmFlag=strcmpi(cfg.phaseTrnsfrm,'yes');


    
%Number of data points per calculation and default number of bins: 
switch method
    case 1 %trial
        NpointsPcalc = cfg.Ncut;
    case 2 %trialTime
        NpointsPcalc = cfg.Nwin;
    case 3 %ensemble
         NpointsPcalc = cfg.Nwin*cfg.Ntr;
end
defaultNbins = max( 5 , min(2*floor(sqrt(NpointsPcalc)/2)+1,33) );

if isfield(cfg,'Nbins')
    varName='Nbins';
    testNBINS = {@(Nbins)isnumeric(Nbins),...
                 @(Nbins)isreal(Nbins),...
                 @(Nbins)isscalar(Nbins),... 
                 @(Nbins)( Nbins>4),...
                 @(Nbins)( Nbins<NpointsPcalc ),...
                 @(Nbins)round(Nbins)==Nbins,...
                 @(Nbins)mod(Nbins,2)==1};
    param={{},{},{},{},{},{},{}};
    mode=['e','e','e','e','e','e','w'];
    execfun={{},{},{},{},{},{},@(Nbins)Nbins-1};
    default=defaultNbins;
    [RESULT, cfg.Nbins] = DPvalidateData(cfg.Nbins,testNBINS,param,mode,execfun,default,varName,funcName);
else
    cfg.Nbins = defaultNbins;
    %cprintf('Magenta','WARNING: Structure field Nbins is missing. Replacing with default value: %d.\n',cfg.Nbins)
    fprintf('WARNING: Structure field Nbins is missing. Replacing with default value: %d.\n',cfg.Nbins)
end 
cfg.binLenPhi = TwoPi/cfg.Nbins; %bin length for phi
cfg.binLenDphi = TwoPi/(cfg.Nbins-1); %bin length for Dphi
cfg.gridPhi= cfg.binLenPhi/2:cfg.binLenPhi:TwoPi-cfg.binLenPhi/2; %[0 2pi] bin grid for phi
cfg.gridPhiEdgs = 0:cfg.binLenPhi:TwoPi;
cfg.gridDphi = - pi:cfg.binLenDphi:pi;%[-pi pi] bin grid for Dphi


if isfield(cfg,'measures')
    varName='measures';
    testMEASURES = { @(measures)iscellstr(measures)||ischar(measures),...
                     @(measures)~isempty(intersect(measures,{'all','IC','PLV','PPC','PLI','SE','CP','MI','PR'})) };
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
% %             structure as constructed here:
MeasNames={'PCI','NCI','ACI','ICI12','ICI21','','';...
    'PLV','','','','','','';...
    'PPC','','','','','','';...
    'PLI','','','','','','';...
    'SE','','','','','','';...
    'CP','','','','','','';...
    'MI','','','','','','';...
    'PRPL','PRcn1','PRcn2','PRcn','PRcd1','PRcd2','PRcd';...
    'CCR','','','','','','';...
    };
% 
% % -NmeasPmeas: the vector of the numbers of submeasures per measure 
% %              in the results' structure as constructed here:
NmeasPmeas = [5 1 1 1 1 1 1 7 1]; %Number of 'measures per measures'
% categories


%Unpack the measures that will be calculated
if strcmpi(cfg.measures,'all')
    measures = ones(1,9);
else
    measures = [any(strcmpi(cfg.measures,'IC')),...
                any(strcmpi(cfg.measures,'PLV')),...
                any(strcmpi(cfg.measures,'PPC')),...
                any(strcmpi(cfg.measures,'PLI')),...
                any(strcmpi(cfg.measures,'SE')),...
                any(strcmpi(cfg.measures,'CP')),...
                any(strcmpi(cfg.measures,'MI')),...
                any(strcmpi(cfg.measures,'PR')),...
                any(strcmpi(cfg.measures,'CCR'))];
end
if measures(8) %PR
    %Check whether PR measures can be calculated: they need at least 1
    %period time windows. Shorter windows or pointwise calculation are not
    %possible
    if TPR < cfg.TcMax
        %cprintf('Magenta','WARNING: PR measures cannot be calculated for a time window shorter than 1 period. PR calculation is excluded.\n')
        fprintf('WARNING: PR measures cannot be calculated for a time window shorter than 1 period. PR calculation is excluded.\n')
        measures(8)=0;
    end
    
end
measInds = find(measures); %Indexes of measures to be calculated
Nmeasures = length(measInds); %number of measures to calculate
if (Nmeasures==0)
    error('Exiting because there are no measures to calculate!')
end


%Measures' parameters structures:
if measures(1)
    if isfield(cfg,'IC')
        varName='IC';
        testIC = {@(IC)isstruct(IC) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.IC] = DPvalidateData(cfg.IC,testIC,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.IC,'Dphi0')
            varName='IC.Dphi0';
            testICDPHI0 = {@(Dphi0)isnumeric(Dphi0),...
                           @(Dphi0)isreal(Dphi0),...
                           @(Dphi0)( Dphi0>0 )&&( Dphi0<=pi/2 ),...
                           @(Dphi0)isscalar(Dphi0) };
            param={{},{},{},{}};
            mode=['e','e','e','e'];
            execfun={{},{},{},{}};
            default=pi/4;
            [RESULT, cfg.IC.Dphi0] = DPvalidateData(cfg.IC.Dphi0,testICDPHI0,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.IC.Dphi0 is missing. Replacing with default value pi/4.\n')
            fprintf('WARNING: Structure field cfg.IC.Dphi0 is missing. Replacing with default value pi/4.\n')
            cfg.IC.Dphi0 = pi/4;
        end
    else
        %cprintf('Magenta','WARNING: Structure field cfg.IC is missing. Replacing with default value cfg.IC.Dphi0=pi/4.\n')
        fprintf('WARNING: Structure field cfg.IC is missing. Replacing with default value cfg.IC.Dphi0=pi/4.\n')
        cfg.IC.Dphi0=pi/4;
    end
end


if measures(8)
    if isfield(cfg,'PR')
        varName='PR';
        testPR = {@(PR)isstruct(PR) };
        param = {{}};
        mode=['e'];
        execfun={{}};
        default=nan;
        [RESULT, cfg.PR] = DPvalidateData(cfg.PR,testPR,param,mode,execfun,default,varName,funcName);
        
        if isfield(cfg.PR,'method')
            varName='PR.method';
            testPRMETHOD = { @(method)ischar(method),...
                             @(method)isvector(method),...
                             @(method)strcmpi(method,'fourier') || strcmpi(method,'iter') };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='fourier';
            [RESULT, cfg.PR.method] = DPvalidateData(cfg.PR.method,testPRMETHOD,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.PR.method is missing. Replacing with default value ''fourier''.\n')
            fprintf('WARNING: Structure field cfg.PR.method is missing. Replacing with default value ''fourier''.\n')
            cfg.PR.method = 'fourier';
        end

        
        if isfield(cfg.PR,'Ngrid')
            varName='PR.Ngrid';
            testPRNGRID = {@(Ngrid)isnumeric(Ngrid),...
                           @(Ngrid)isreal(Ngrid),...
                           @(Ngrid)( Ngrid>0 ),...
                           @(Ngrid)round(Ngrid)==Ngrid,...
                           @(Ngrid)isscalar(Ngrid) };
            param={{},{},{},{},{}};
            mode=['e','e','e','e','e'];
            execfun={{},{},{},{},{}};
            default=100;
            [RESULT, cfg.PR.Ngrid] = DPvalidateData(cfg.PR.Ngrid,testPRNGRID,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.PR.Ngrid is missing. Replacing with default value 100.\n')
            fprintf('WARNING: Structure field cfg.PR.Ngrid is missing. Replacing with default value 100.\n')
            cfg.PR.Ngrid = 100;
        end
        
        if isfield(cfg.PR,'or')
            varName='PR.or';
            testPROR = {@(or)isnumeric(or),...
                           @(or)isreal(or),...
                           @(or)( or>0 ),...
                           @(or)round(or)==or,...
                           @(or)isscalar(or) };
            param={{},{},{},{},{}};
            mode=['e','e','e','e','e'];
            execfun={{},{},{},{},{}};
            default=10;
            [RESULT, cfg.PR.or] = DPvalidateData(cfg.PR.or,testPROR,param,mode,execfun,default,varName,funcName);
        else
            %cprintf('Magenta','WARNING: Structure field cfg.PR.or is missing. Replacing with default value 10.\n')
            fprintf('WARNING: Structure field cfg.PR.or is missing. Replacing with default value 10.\n')
            cfg.PR.or = 10;
        end
        
    else
        %cprintf('Magenta','WARNING: Structure field cfg.PR is missing. Replacing with default values cfg.PR.method=''fourier'', cfg.PR.Ngrid=100, cfg.PR.or=10.\n')
        fprintf('WARNING: Structure field cfg.PR is missing. Replacing with default values cfg.PR.method=''fourier'', cfg.PR.Ngrid=100, cfg.PR.or=10.\n')        
        cfg.PR.method='fourier';
        cfg.PR.Ngrid=100;
        cfg.PR.or=10;
    end
end



        
%Calculate command strings
if (method==1) %trial
    
    funCommands  = {'[PC.(''PCI'')(iTr), PC.(''NCI'')(iTr), PC.(''ACI'')(iTr), PC.(''ICI12'')(iTr), PC.(''ICI21'')(iTr)] = DPcalcIC(thisDphi,cfg.TcMax,cfg.IC.Dphi0);',...
                    'PC.(''PLV'')(iTr) = DPcalcPLV(thisDphi);', ...
                    'PC.(''PPC'')(iTr) = DPcalcPPC(thisDphi);', ...
                    'PC.(''PLI'')(iTr) = DPcalcPLI(thisDphi);', ...
                    'PC.(''SE'')(iTr) = DPcalcSE(thisDphi,cfg.Nbins,cfg.gridDphi);', ...                   
                    'PC.(''CP'')(iTr) = DPcalcCP(thisPhi,cfg.Nbins,cfg.gridPhiEdgs);', ...
                    'PC.(''MI'')(iTr) = DPcalcMIphase(thisPhi,cfg.Nbins,cfg.gridPhi);', ...
                    '[PC.(''PRPL'')(iTr),PC.(''PRcn1'')(iTr), PC.(''PRcn2'')(iTr), PC.(''PRcn'')(iTr),PC.(''PRcd1'')(iTr), PC.(''PRcd2'')(iTr), PC.(''PRcd'')(iTr) ]= DPcalcPR(thisPhi,cfg.fs,cfg.phaseTrnsfrmFlag,cfg.PR.method,cfg.PR.Ngrid,cfg.PR.or);', ...
                    'PC.(''CCR'')(iTr) = DPcalcCCorr(thisPhi);', ... 
                    };

elseif (method==2) %trialTime 
 
     funCommands  = {'[PC.(''PCI'')(iT,iTr), PC.(''NCI'')(iT,iTr), PC.(''ACI'')(iT,iTr), PC.(''ICI12'')(iT,iTr), PC.(''ICI21'')(iT,iTr)] = DPcalcIC(thisDphi,cfg.TcMax,cfg.IC.Dphi0);',...
                     'PC.(''PLV'')(iT,iTr) = DPcalcPLV(thisDphi);', ...
                     'PC.(''PPC'')(iT,iTr) = DPcalcPPC(thisDphi);', ...
                     'PC.(''PLI'')(iT,iTr) = DPcalcPLI(thisDphi);', ...
                     'PC.(''SE'')(iT,iTr) = DPcalcSE(thisDphi,cfg.Nbins,cfg.gridDphi);', ...
                     'PC.(''CP'')(iT,iTr) = DPcalcCP(thisPhi,cfg.Nbins,cfg.gridPhiEdgs);', ...
                     'PC.(''MI'')(iT,iTr) = DPcalcMIphase(thisPhi,cfg.Nbins,cfg.gridPhi);', ...
                     '[PC.(''PRPL'')(iT,iTr),PC.(''PRcn1'')(iT,iTr), PC.(''PRcn2'')(iT,iTr), PC.(''PRcn'')(iT,iTr),PC.(''PRcd1'')(iT,iTr), PC.(''PRcd2'')(iT,iTr), PC.(''PRcd'')(iT,iTr) ] = DPcalcPR(thisPhi,cfg.fs,cfg.phaseTrnsfrmFlag,cfg.PR.method,cfg.PR.Ngrid,cfg.PR.or);',...
                     'PC.(''CCR'')(iT,iTr) = DPcalcCCorr(thisPhi);', ...
                     };

                
else %ensemble 
    
    funCommands  = {'[PC.(''PCI'')(iT), PC.(''NCI'')(iT), PC.(''ACI'')(iT), PC.(''ICI12'')(iT), PC.(''ICI21'')(iT)] = DPcalcICensemble(thisDphi,cfg.NtcMax,cfg.IC.Dphi0,thisNwin,cfg.Ntr);',...
                    'PC.(''PLV'')(iT) = DPcalcPLV(thisDphi);', ...
                    'PC.(''PPC'')(iT) = DPcalcPPC(thisDphi);', ...
                    'PC.(''PLI'')(iT) = DPcalcPLI(thisDphi);', ...
                    'PC.(''SE'')(iT) = DPcalcSE(thisDphi,cfg.Nbins,cfg.gridDphi);', ...                                       
                    'PC.(''CP'')(iT) = DPcalcCP(thisPhi,cfg.Nbins,cfg.gridPhiEdgs);', ...
                    'PC.(''MI'')(iT) = DPcalcMIphase(thisPhi,cfg.Nbins,cfg.gridPhi);',...
                    '[PC.(''PRPL'')(iT),PC.(''PRcn1'')(iT), PC.(''PRcn2'')(iT), PC.(''PRcn'')(iT),PC.(''PRcd1'')(iT), PC.(''PRcd2'')(iT), PC.(''PRcd'')(iT) ]=  DPcalcPRensemble(thisPhi,cfg.fs,cfg.phaseTrnsfrmFlag,cfg.PR.method,cfg.PR.Ngrid,cfg.PR.or,thisNwin,cfg.Ntr);'...
                    'PC.(''CCR'')(iT) = DPcalcCCorr(thisPhi);', ... 
                    };
end      

%Select the function commands of the measures to be calculated
thisfunCommands = funCommands(measInds); 

%If a transformation from proto- to true phases is need...
if cfg.phaseTrnsfrmFlag
    %...append it at the beginning of the commands to be executed
    thisfunCommands = ['thisPhi(:,1) = co_fbtransf1(thisPhi(:,1));  thisPhi(:,2) = co_fbtransf1(thisPhi(:,2));',thisfunCommands];
end
Ncomnds = length(thisfunCommands);

%Hilbert transform and optional cutting of edges if input data is a raw
%signal
if strcmpi(cfg.domain,'real')
    phi=zeros(cfg.Ncut,2,cfg.Ntr);
    for iT=1:Ntr; 
        [phi(:,:,iT), dummy, dummy2] = DPco_hilbproto(x(:,:,iT),cfg.cutTailsPoints);
    end
else
    phi=x;
    clear x;
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
        if strcmpi(cfg.stats.corrMultComp,'BONF')
            switch method
                case 1
                    defaultNperm = cfg.Ntr*defaultNperm;
                case 2
                    defaultNperm =  cfg.Ntr*cfg.Ncalc*defaultNperm;
                case 3
                    defaultNperm =  cfg.Ncalc*defaultNperm;
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
