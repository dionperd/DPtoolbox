function [C, cfg, statsRes] = DPstateCoupling(cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas, x,y,z)

%%This function prepares data and configuration struture for DPstateCoupling
%

%Inputs: 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-cfg: configuration structure with fields
%   -fs: sampling frequency in herz, positive real scalar
%   -fc: central frequency of the signals in herz, 2 (or 3) element vector of
%        positive real values
%   -time: time vector in secs, real valued vector, default=[0:N-1]/fs
%   -domain: 'real' or 'state', string, if domain=='real'
%   -Dx: number of time delay dimensions of x,
%   -Dy: number of time delay dimensions of y,
%   -Dz (optionally): number of time delay dimensions of z,
%   -tauX: a vector of Dx-1 embedding time delays of x, in samples
%   -tauY: a vector of Dy-1 embedding time delays of y, in samples
%   -tauZ (optionally): a vector of Dz-1 embedding time delays of z, in samples
%   -u: range of time lags to be tested, in samples, vector of real numbers >=0
%       the test is calculated symmetrically for time lages +/- u(i), 
%       default, depending on the measure: 
%        -0 for MI, TI, NI, CCD, 
%        -and 1/fs, i.e. one sampling point, for TE, PTE, MP 
%   -Wth: the Theiler window to be used for the nearest neighbors
%         statistics in samples
%--------------------------------------------------------------------------
%Optionally if domain=='real' then a time delay embedding is needed:
%        Then, Np = min([Nx,Ny,Nz]);
%--------------------------------------------------------------------------
%   -method: one of 
%            1. 'trial', for estimation per trial
%            2. 'trialTime', for estimation per trial with a sliding time 
%               window
%            3. 'ensemble', for estimation across trials, either pointwise
%            or with a time window
%--------------------------------------------------------------------------
%   -zscoreTrnsfrm: 'yes' or 'no',string, default='yes' 
%                but only if domain='real' 
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
%%   -measures: the coupling measures that should be calculated 
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
%  -X: embedded pointset of signal x, matrix (Np,Dx,Ntr)
%  -Y: embedded pointset of signal y, matrix (Np,Dy,Ntr)
%  Optionally:
%  -Z: embedded pointset of signal z, matrix (Np,Dz,Ntr)


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


tic

%Main calculation
switch method
    
    
    case 1 %'trial'
        
        %The dimensions of the output of each measure
        cfg.Nout(1) = cfg.Nu;
        cfg.Nout(2) = cfg.Ntr;
        
        %The time vector of the output
        cfg.timeOut = cfg.timeP;
        
        %Initialize measure vectors:
        %For each measure we will calculate...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Nu,cfg.Ntr);
                if Ntr>1
                    C.trialMean.(MeasNames{measInds(iM),iSubM})=zeros(1,cfg.Nu);
                end
            end
        end
        
        
        %The calculation loop:
        
        %For each time lag u...
        for iU=1:cfg.Nu;
            
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
                        C.trialMean.(MeasNames{measInds(iM),iSubM})(iU) = nanmean( C.(MeasNames{measInds(iM),iSubM})(iU,:) );
                    end
                end
            end
            
        end
        
    case 2 %'trialTime'
        
        %The dimensions of the output of each measure
        cfg.Nout(1) = cfg.Ncalc;
        cfg.Nout(2) = cfg.Nu;
        cfg.Nout(3) = cfg.Ntr;
        
        %The time vector of the output
        cfg.timeOut = cfg.timeCalc;
        
        %Initialize measure vectors:
        %For each measure we will calculate...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,cfg.Nu,cfg.Ntr);
            end
        end
        
        

        %The calculation loop:
        
        %Initialize time window centers and edges for this u
        cfg.timePointWinEdgs = nan(cfg.Ncalc,2);
        cfg.timeWinEdgs = nan(cfg.Ncalc,2);
        
        %For each time lag u...
        for iU=1:cfg.Nu;

            %For each calculation/window step...
            for iT=1:cfg.Ncalc;
                
                if (iU==1)
                    %...calculate indexes for:
                    cfg.timePointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
                    cfg.timePointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.Np);%ending window point
                    cfg.timeWinEdgs(iT,:) = cfg.timeP(cfg.timePointWinEdgs(iT,:)); %same points in time...
                end
                
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
                    
                    clear thisPhi thisDphi;
                    
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
                        C.trialMean.(MeasNames{measInds(iM),iSubM})(:,iU) = nanmean( squeeze(C.(MeasNames{measInds(iM),iSubM})(:,iU,:)), 2 );
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
                        %...and for each time lag u...
                        for iU=1:cfg.Nu;
                            %...upsample from Ncalc points to Np:
                            C.trialMean.(MeasNames{measInds(iM),iSubM})(:,iU) = resample( oldC.(MeasNames{measInds(iM),iSubM})(:,iU), cfg.Np, cfg.Ncalc );
                        end
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
        cfg.Nout(2) = cfg.Nu;
        
        %The time vector of the output
        cfg.timeOut = cfg.timeCalc;
        
        %Initialize measure vectors:
        %For each measure we will calculate...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                %...initialize with zeros
                C.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,cfg.Nu);
            end
        end
        
        %Initialize time window centers and edges
        cfg.timePointWinEdgs = nan(cfg.Ncalc,2);
        cfg.timeWinEdgs = nan(cfg.Ncalc,2);
        
        
        %For each time lag u...
        for iU=1:cfg.Nu;
            
            %For each calculation/window step...
            for iT=1:cfg.Ncalc;
                
                if (iU==1)
                    %...calculate indexes for:
                    cfg.timePointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
                    cfg.timePointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.Np);%ending window point
                    cfg.timeWinEdgs(iT,:) = cfg.timeP(cfg.timePointWinEdgs(iT,:)); %same points in time...
                end

                
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
            
        end
        
        %Interpolate through upsampling if required:
        if strcmpi(cfg.upsample,'yes')
            oldC = C; %temporary copy
            %For each measure we calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure ...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...and for each time lag u...
                    for iU=1:cfg.Nu;
                        %...upsample from Ncalc points to Np:
                        C.(MeasNames{measInds(iM),iSubM}) = resample( oldC.(MeasNames{measInds(iM),iSubM}), cfg.Np, cfg.Ncalc );
                    end
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
        dataSurr = cfg.stats.surrfun([x y],2,cfg.Ntr);
        xSurr = dataSurr(:,1);
        ySurr = dataSurr(:,2);
        clear dataSurr;
        %...and calculate the measures with exactly the same setting
        %except for the surrogate statistics...
        disp(['...calculating measures for surrogate data set ',num2str(iST),'/',num2str(cfg.stats.Nperm),'...'])
        [Csurr{iST}, ~, ~] = DPstateCoupling(cfgSurr, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas,xSurr,ySurr,z); 
        
    end
    
    
    %Initialize statistics results structure
    statsRes = struct();
    
    %...initialize measure vectors:
    %For each measure we have calculated...
    for iM = 1:Nmeasures;
        %...and for each of the sumbmeasures of this measure...
        for iSubM = 1:NmeasPmeas(measInds(iM));
            
            disp(['...calculating statistics for measure ',MeasNames{measInds(iM),iSubM},'...'])
            
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
 


