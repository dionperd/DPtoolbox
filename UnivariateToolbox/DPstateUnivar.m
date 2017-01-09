function [C, cfg, statsRes] = DPstateUnivar(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds,x,y,z)

%%This function prepares data and configuration struture for DPstateCoupling
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
%--------------------------------------------------------------------------
%  -method: 1,2,or 3 depending on whether cfg.method =  'trial', 
%           'trialTime', or 'ensemble' respectively 
%  -measInds: the indexes of the measures to be calculated
%  -Nmeasures: how many measures are to be calculated
%  -thisfunCommands: cell of the command strings for the calculation of
%                    the selected measures
%  -X: embedded pointset of signal x, matrix (Np,D,Ntr)



%  Outputs:

%  -C: results' structure with fields with the names of the calculated measures
%       of size depending on the method:
%       1. vector 1 x M for 'method='trial'
%       2. matrix cfg.Ncalc x M for 'method='trialTime'
%       3. vector Nout x 1  for method='ensemble'
% if method == 'trial' or 'trialTime', then mean values across trials are
% given in a substructure C.trialMean as:
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




%Define a constant:
TwoPi = 2*pi;

%Cell of names of result's structure
MeasNames={'CD';...
            };
NmeasPmeas = [1]; %Number of 'measures per measures' categories
NoutsPmeas = [1;...
                ];


tic

%Main calculation
switch method
    
    
    case 1 %'trial'
        
        %The dimensions of the output of each measure
        cfg.Nout(1) = cfg.Ntr;
        
        %The time vector of the output
        cfg.timeOut = cfg.timeP;
        
        %Initialize measure vectors:
        %For each measure we will calculate...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                Nouts = cfg.NoutsPmeas(measInds(iM),iSubM);
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(Nouts,cfg.Ntr);
                if Ntr>1
                    C.trialMean.(MeasNames{measInds(iM),iSubM})=zeros(1,Nouts);
                end
            end
        end
        
        
        %The calculation loop:
        
        %...for each trial...
        for iTr=1:cfg.Ntr;
            
            %disp(['...calculating trial',num2str(iT),'/',num2str(cfg.Ntr),'...'])
            
            %...calculate the selected measures for this trial
            for iC = 1:Ncomnds;
                try
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
                    C.trialMean.(MeasNames{measInds(iM),iSubM}) = nanmean( C.(MeasNames{measInds(iM),iSubM}) );
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
                Nouts = cfg.NoutsPmeas(measInds(iM),iSubM);
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,Nouts,cfg.Ntr);
            end
        end
        
        
        
        %The calculation loop:
        
        %Initialize time window centers and edges for this u
        cfg.timePointWinEdgs = nan(cfg.Ncalc,2);
        cfg.timeWinEdgs = nan(cfg.Ncalc,2);
        
        
        %For each calculation/window step...
        for iT=1:cfg.Ncalc;
            
            %...calculate indexes for:
            cfg.timePointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
            cfg.timePointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.Np);%ending window point
            cfg.timeWinEdgs(iT,:) = cfg.timeP(cfg.timePointWinEdgs(iT,:)); %same points in time...
            
            
            %...for each trial...
            for iTr=1:cfg.Ntr;
                
                %                     %...get the time series of this trial and time window
                %                     thisX = squeeze( x( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),:,iTr ) );
                %                     thisY = squeeze( y( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),:,iTr ) );
                
                %...calculate the selected measures for this trial and time point
                for iC = 1:Ncomnds;
                    try
                        eval(thisfunCommands{iC});
                    catch
                        keyboard
                    end
                end
                
                
            end
            
            
            %Calculate the mean of all trials:
            if cfg.Ntr>1
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...calculate the mean value (excluding possible nan
                        %values)
                        C.trialMean.(MeasNames{measInds(iM),iSubM}) = nanmean( squeeze(C.(MeasNames{measInds(iM),iSubM})), 2 );
                    end
                end
            end
            
        end
            
            %Interpolate through upsampling if required:
            if strcmpi(cfg.upsample,'yes')
                oldC = C.trialMean; %temporary copy
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...upsample from Ncalc points to Np:
                        C.trialMean.(MeasNames{measInds(iM),iSubM}) = resample( oldC.(MeasNames{measInds(iM),iSubM}), cfg.Np, cfg.Ncalc ); 
                    end
                end
                cfg.Nout(1) = cfg.Np;
                cfg.timeOut = cfg.timeP;
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
                Nouts = cfg.NoutsPmeas(measInds(iM),iSubM);
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,Nouts);
            end
        end
        
        %Initialize time window centers and edges
        cfg.timePointWinEdgs = nan(cfg.Ncalc,2);
        cfg.timeWinEdgs = nan(cfg.Ncalc,2);
        
        
        %For each calculation/window step...
        for iT=1:cfg.Ncalc;
            
            %...calculate indexes for:
            cfg.timePointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
            cfg.timePointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.Np);%ending window point
            cfg.timeWinEdgs(iT,:) = cfg.timeP(cfg.timePointWinEdgs(iT,:)); %same points in time...
            
            
            %                     %...get the time series of this trial and time window
            %                     thisX = squeeze( x( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),:,: ) );
            %                     thisY = squeeze( y( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),:,: ) );
            
            %...calculate the selected measures for this trial and time point
            for iC = 1:Ncomnds;
                try
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
                    %...upsample from Ncalc points to Np:
                    C.(MeasNames{measInds(iM),iSubM}) = resample( oldC.(MeasNames{measInds(iM),iSubM}), cfg.Np, cfg.Ncalc ); 
                end
            end
            cfg.Nout(1) = cfg.Np;
            cfg.timeOut = cfg.timeP;
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
        
        toc
        
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
    
    %...for each pair of the surrogate time series...
    for iST=1:cfg.stats.Nperm;
        
        %...create it...
        xSurr = cfg.stats.surrfun(x,1,cfg.Ntr);

        %...and calculate the measures with exactly the same setting
        %except for the surrogate statistics...
        disp(['...calculating measures for surrogate data set ',num2str(iST),'/',num2str(cfg.stats.Nperm),'...'])
        [Csurr{iST}, ~, ~] = DPstateUnivar(cfgSurr, method, measInds, Nmeasures, thisfunCommands, Ncomnds,xSurr); 
        
    end
    
    
    %Initialize statistics results structure
    statsRes = struct();
    
    %...initialize measure vectors:
    %For each measure we have calculated...
    for iM = 1:Nmeasures;
        %...and for each of the sumbmeasures of this measure...
        for iSubM = 1:NmeasPmeas(measInds(iM));
                            
            disp(['...calculating statistics for measure ',MeasNames{measInds(iM),iSubM},'...'])
            
            Nouts = cfg.NoutsPmeas(measInds(iM),iSubM);
            
            %...select method...
            switch method
                
                case 1 %trial
                    
                    
                    %...create a temporary matrix of the measure of the surrogate dataset...
                    thisCsurr =  zeros(cfg.Ntr,cfg.stats.Nperm);
                    for iST=1:cfg.stats.Nperm;
                        %...for each permutation...
                        thisCsurr(:,iST) = Csurr{iST}.(MeasNames{measInds(iM),iSubM});
                    end
                    
                    %...calculate the specific statistic for each trial...
                    [pointStat,  ~] = DPcalcStat(C.(MeasNames{measInds(iM),iSubM}),thisCsurr,pointStatMethod, multiStatfun);
                    clear thisPCsurr;
                    
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
                        thisPCmeanSurr =  zeros(1,cfg.stats.Nperm);
                        for iST=1:cfg.stats.Nperm;
                            thisCmeanSurr(iST) = Csurr{iST}.trialMean.(MeasNames{measInds(iM),iSubM});
                        end
                        
                        %...calculate the specific statistic...
                        [pointStat,  ~] = DPcalcStat(C.trialMean.(MeasNames{measInds(iM),iSubM}),thisCmeanSurr,pointStatMethod, multiStatfun);
                        clear thisPCmeanSurr;
                        
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
                    thisCsurr =  zeros(cfg.Nout(1),cfg.Ntr,cfg.stats.Nperm);
                    for iST=1:cfg.stats.Nperm;
                        %...for each trial...
                        thisCsurr(:,:,iST) = Csurr{iST}.(MeasNames{measInds(iM),iSubM});
                    end
                    
                    %...calculate the specific statistic for each trial...
                    [pointStat,  multiStat] = DPcalcStat(C.(MeasNames{measInds(iM),iSubM}),thisCsurr,pointStatMethod, multiStatfun);
                    clear thisPCsurr;
                    
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
                        thisCmeanSurr =  zeros(cfg.Nout(1),cfg.stats.Nperm);
                        for iST=1:cfg.stats.Nperm;
                            thisCmeanSurr(:,iST) = Csurr{iST}.trialMean.(MeasNames{measInds(iM),iSubM});
                        end
                        
                        %...calculate the specific statistic...
                        [pointStat,  multiStat] = DPcalcStat(C.trialMean.(MeasNames{measInds(iM),iSubM}),thisCmeanSurr,pointStatMethod, multiStatfun);
                        clear thisPCmeanSurr;
                        
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
                    thisCsurr =  zeros(cfg.Nout(1),cfg.stats.Nperm);
                    for iST=1:cfg.stats.Nperm;
                        thisCsurr(:,iST) = Csurr{iST}.(MeasNames{measInds(iM),iSubM});
                    end
                    
                    %...calculate the specific statistic...
                    [pointStat,  multiStat] = ...
                        DPcalcStat(C.(MeasNames{measInds(iM),iSubM}),thisCsurr,pointStatMethod, multiStatfun);
                    clear thisPCsurr;
                    
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
    
end

disp('DONE!')



function [pointStat,  multiStat] = DPcalcStat(mVal,surrVal,pointStatMethod,multiStatfun)

%This function calculates a point statistic according to pointStatMethod
%and a multivariate one, according to the function handle multiVarStat


SIZEm = size(mVal); %the size of the measurement values
N = numel(mVal); %total number of the measurement values
SIZEsurr = size(surrVal); %the size of the surrogate values
Nsurr = SIZEsurr(end); %the number of surrogate values per data point

% %Initialize
% pointStat.mval = zeros(N,1);
% pointStat.surrVal = zeros(N,Nsurr);

%Vectorize the statistic
mVal = mVal(:);

%Collapse extra dimensions of surrVal, if any
if length(SIZEsurr)>2
    surrVal  = reshape(surrVal,[N, Nsurr]);
end

%Calculate the mean...
pointStat.mean = reshape(mean(surrVal,2), SIZEm );
%...the standard deviation...
pointStat.std = reshape(std(surrVal,0,2), SIZEm );
%...the most common value...
pointStat.mode = reshape(mode(surrVal,2), SIZEm );
%...the median...
pointStat.median = reshape(median(surrVal,2), SIZEm );
%...the min...
pointStat.min = reshape(min(surrVal,[],2), SIZEm );
%...and the max...
pointStat.max = reshape(max(surrVal,[],2), SIZEm );
%...of the surrogate values

if strcmpi(pointStatMethod,'z')
    
    pointStat.mVal = reshape( ( mVal -pointStat.mean ) ./pointStat.std, SIZEm );
    
    pointStat.surrVal = reshape( ( surrVal - repmat(pointStat.mean,[1 Nsurr]) ) ./repmat(pointStat.std,[1,Nsurr]), SIZEsurr );
    
elseif strcmpi(pointStatMethod,'t')
    
    pointStat.mVal = reshape( sqrt(Nsurr) * ( mVal -pointStat.mean ) ./pointStat.std, SIZEm );
    
    pointStat.surrVal = reshape( sqrt(Nsurr) * ( surrVal - repmat(pointStat.mean,[1 Nsurr]) ) ./repmat(pointStat.std,[1,Nsurr]), SIZEsurr );
    
else
    pointStat.mVal = reshape( mVal, SIZEm );
    
    pointStat.surrVal = reshape( surrVal, SIZEsurr );
end

mVal = reshape( mVal, SIZEm );
surrVal = reshape( surrVal, SIZEsurr );

if ~isempty(multiStatfun) && (SIZEm(1)>1)   
    multiStat = multiStatfun.fun(mVal,surrVal,pointStat,multiStatfun.params);  
else
    multiStat=[];
end



function [p, pUnCorr, sign] = DPcalcSignif(stat,alpha,tail,corrMultComp)

%This function decides about the significance of a specific statistic in a
%surrogate test

%Inputs:
%   -stat: a structure containing the statistic, it contains the fields:
%          -mVal: the measurement value array
%          -surrVal: the surrogate values' array
%   -alpha: the alpha value
%   -tail: 1 or 2, for one or two tailed test
%   -corrMultComp: method of correction for multiple comparisons, '','BONF' or
%                  'FDR', default:'' for no correction

%Outputs: 
%   -p: an array of p-values after correction for multiple comparisons
%   -pUnCorr: an array of p-values before correction for multiple
%             comparisons (p=pUnCorr if there is no correction)
%   -sign: an array of 1s and 0s for the significant and non significant
%          data points, respectively

%p, pUnCorr and sign have the same structure as mVal.



SIZEm = size(stat.mVal); %the size of the measurement values
N = numel(stat.mVal); %total number of the measurement values
SIZEsurr = size(stat.surrVal); %the size of the surrogate values
Nsurr = SIZEsurr(end); %the number of surrogate values per data point
NsurrS = repmat(Nsurr,[N,1]); %minimum p value

%Initialize outputs:
p=nan(N,1);
pUnCorr=nan(N,1);
sign=nan(N,1);

%Vectorize the statistic
stat.mVal = stat.mVal(:);
mVals = repmat(stat.mVal, [1, Nsurr]);

%Collapse extra dimensions of surrVal, if any
if length(SIZEsurr)>2
    stat.surrVal  = reshape(stat.surrVal,[N, Nsurr]);
end

if (tail==1) %For one tailed test...
    
    %...calculate the probability of each data point to be smaller than the
    %surrogate values, with a minimum of 1/Nsurr...
    pUnCorr = max( 1, sum( mVals<stat.surrVal ,2) )/ Nsurr;
    
    
else %For two tailed test...
    
    %...calculate the probability of each data point to be larger or smaller than the
    %surrogate values, with a minimum of 1/Nsurr...
    pUnCorr = max( 1, min( sum( mVals>stat.surrVal ,2),  sum( mVals<stat.surrVal ,2) ) )/Nsurr;
    
    %...divide alpha value by 2 for two-tailed test...
    alpha = alpha/2;
    
end


if strcmpi(corrMultComp,'BONF') %BONF
    
    %...p equals pUnCorr multiplied by number of comparisons...
    p = min(pUnCorr*N, 1);
    
    
elseif strcmpi(corrMultComp,'FDR') % 'FDR'
    
    %...sort p values...
    [pSort, sortIndx] = sort(pUnCorr);
    
    %...FDR correction...
    NN=1./(1:N); %a useful constant
    p= sum(NN)*pSort.*(N.*NN); %correction
    p = min(p, 1);
    clear pSort;
    
    %...unsort back...
    [~,unsortIndx]=sort(sortIndx);
    p = p(unsortIndx);
    
    
else %...if there is no correction for multiple comparisons...
    
    %...p equals pUnCorr...
    p=pUnCorr;
    
    
end
    
%Calculate significance...
sign = p<=alpha;

%Reshape results back to mVal size:
p=reshape(p,SIZEm);
pUnCorr=reshape(pUnCorr,SIZEm);
sign=reshape(sign,SIZEm);



    


