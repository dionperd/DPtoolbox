function [PC, Dphi, cfg, statsRes] = DPphaseCoupling(phi0, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas)

%This function calculates several phase coupling measures optionaly in a
%time resolved manner and/or accross trials


%  Inputs:

%  -phi0: phases of the signals
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


%for cfg:
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

%  -PC: results cell for the K different nm calculations
%       of structures with fields with the names of the calculated measures
%       of size depending on the method:
%       1. vector 1 x M for 'method='trial'
%       2. matrix cfg.Ncalc x M for 'method='trialTime'
%       3. vector Ncut x Nout  for method='ensemble'
% if method == 'trial' or 'trialTime', then mean values across trials are
% given in a substructure PC.trialMean as:
%       1. scalar for 'method='trial'
%       2. vector cfg.Nout x 1 for 'trialTime'
%  -Dphi: cell of the phase difference in the interval [-pi pi], of
%         dimensions Ncut x M, for the K different nm calculations
%  -cfg: configuration structure corrected and complemented
%  -statsRes: results cell structures for the surrogate statistics' test,
%             same structure as PC, with fieds
%             -poinStat: values of the statistic per data point
%             -multiStat: values of multivariate statistic
%             -p: p values
%             -pUncorr: p values without correction for multiple
%                       comparisons
%             -sign: flags of value 1 or 0 for significant or non
%                    significant results
%



%Define a constant:
TwoPi = 2*pi;


[K dummy]=size(cfg.nm); %number of different nm calculations

%Initialize result cells:
DphiCell = cell(max(cfg.nm(:,1)),max(cfg.nm(:,2)));
PCcell = cell(max(cfg.nm(:,1)),max(cfg.nm(:,2)));


tic
%Main loop for each n,m combination
for iNM=1:K;
    
    %This nm:
    nm=cfg.nm(iNM,:);
    
    disp(['Calculating measures for ',num2str(nm(1)),':',num2str(nm(2)),' coupling...'])
    
    %Calculate phi for n:m combination and the respective Dphi
    phi=phi0;
    for iP=1:2;
        if nm(iP)>1
            phi(:,iP,:) = nm(iP)*phi0(:,iP,:);
        end
    end
    Dphi = squeeze( mod( phi(:,1,:) - phi(:,2,:) , TwoPi ) );
    Dphi(Dphi>pi) = Dphi(Dphi>pi) - TwoPi;
    phi=mod(phi,TwoPi);
    
    %Main calculation
    switch method
        case 1 %'trial'
            
            %The dimensions of the output of each measure
            cfg.Nout(1) = 1;
            cfg.Nout(2) = cfg.Ntr;
            
            %The time vector of the output
            cfg.timeOut = cfg.timeCut;
            
            %Initialize measure vectors:
            %For each measure we will calculate...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...initialize with zeros
                    PC.(MeasNames{measInds(iM),iSubM})=zeros(1,cfg.Ntr);
                end
            end
            
            
            %The calculation loop:
            %For each trial...
            for iTr=1:cfg.Ntr;
                
                %disp(['...calculating trial ',num2str(iT),'/',num2str(cfg.Ntr),'...'])
                
                %...get the phases and phase differences of this trial
                thisPhi = squeeze( phi( :,:,iTr ) );
                thisDphi = Dphi(:,iTr);
                
                %...calculate the selected measures for this trial
                for iC = 1:Ncomnds;
                    try
                        eval(thisfunCommands{iC});
                    catch
                        keyboard;
                    end
                end
                
                clear thisPhi thisDphi;
                
            end
            
            %Calculate the mean of all trials:
            %if cfg.Ntr>1
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...calculate the mean value (excluding possible nan
                        %values)
                        %PC.trialMean.(MeasNames{measInds(iM),iSubM}) = nanmean( PC.(MeasNames{measInds(iM),iSubM}) );
                        PC.trialMean.(MeasNames{measInds(iM),iSubM}) = mean( PC.(MeasNames{measInds(iM),iSubM}) );
                    end
                end
           % end
            
            
        case 2 %'trialTime'
            
            %The dimensions of the output of each measure
            cfg.Nout(1) = cfg.Ncalc;
            cfg.Nout(2) = cfg.Ntr;
                        
            %The time vector of the output
            cfg.timeOut = cfg.timeCalc;
            
            %Initialize measure vectors:
            %For each measure we will calculate...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...initialize with zeros
                    PC.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,cfg.Ntr);
                end
            end
            
            %Initialize time window centers and edges
            cfg.timePointWinEdgs = zeros(cfg.Ncalc,2);
            cfg.timeWinEdgs = zeros(cfg.Ncalc,2);
            
            
            %The calculation loop:
            %For each trial...
            for iTr = 1:cfg.Ntr;
                
                %For each calculation/window step...
                for iT=1:cfg.Ncalc;
                    
                    %...calculate indexes for:
                    cfg.timePointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
                    cfg.timePointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.Ncut);%ending window point
                    cfg.timeWinEdgs(iT,:) = cfg.timeCut(cfg.timePointWinEdgs(iT,:)); %same points in time...
                    
                    %...get the phases and phase differences of this trial and time window
                    thisPhi(:,1) = squeeze( phi( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),1,iTr ) );
                    thisPhi(:,2) = squeeze( phi( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),2,iTr ) );
                    thisDphi     = squeeze( Dphi(cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2), iTr ) );
                    
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
            %if cfg.Ntr>1
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...calculate the mean value (excluding possible nan
                        %values)
                        %PC.trialMean.(MeasNames{measInds(iM),iSubM}) = nanmean( PC.(MeasNames{measInds(iM),iSubM}), 2 );
                        PC.trialMean.(MeasNames{measInds(iM),iSubM}) = mean( PC.(MeasNames{measInds(iM),iSubM}), 2 );
                    end
                end
            %end
            
            %Interpolate through upsampling if required:
            if strcmpi(cfg.upsample,'yes')
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...upsample from Ncalc points to Ncut:
                        PC.trialMean.(MeasNames{measInds(iM),iSubM}) = resample( PC.trialMean.(MeasNames{measInds(iM),iSubM}), cfg.Ncut, cfg.Ncalc );
                    end
                end
                cfg.Nout(1) = cfg.Ncut;
                cfg.timeOut = cfg.timeCut;
            end
            
            %Smooth through convolution with a window function if required:
            if isfield(cfg,'smoothWin')
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...smooth:
                        PC.trialMean.(MeasNames{measInds(iM),iSubM}) = conv( PC.trialMean.(MeasNames{measInds(iM),iSubM}), cfg.smoothWin, 'same' );
                    end
                end
            end
            
            
        otherwise %'ensemble' 
            
            %The dimensions of the output of each measure
            cfg.Nout(1) = cfg.Ncalc;
            cfg.Nout(2) = 1;
            
            %The time vector of the output
            cfg.timeOut = cfg.timeCalc;
            
            %Initialize measure vectors:
            %For each measure we will calculate...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...initialize with zeros
                    PC.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,1);
                end
            end
            
            %Initialize time window centers and edges
            cfg.timePointWinEdgs = zeros(cfg.Ncalc,2);
            cfg.timeWinEdgs = zeros(cfg.Ncalc,2);
            
            
            %The calculation loop:
            %For each calculation/window step...
            for iT=1:cfg.Ncalc;
                
                %...calculate indexes for:
                cfg.timePointWinEdgs(iT,1) = max(cfg.timeCalcIND(iT)-cfg.Nwin2,1); %starting window point
                cfg.timePointWinEdgs(iT,2) = min(cfg.timeCalcIND(iT)+cfg.Nwin2,cfg.Ncut);%ending window point
                cfg.timeWinEdgs(iT,:) = cfg.timeCut(cfg.timePointWinEdgs(iT,:)); %same points in time...
                
                %...get the phases and phase differences of this time
                %window for all trials and vectorize them
                temp  = squeeze( phi( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),1,: ) );
                thisPhi(:,1) = temp(:);
                temp =  squeeze( phi( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),2,: ) );
                thisPhi(:,2) = temp(:);
                temp = squeeze( Dphi( cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2), : ) );
                thisNwin = size(temp,1);
                thisDphi     = temp(:);
                
                %...calculate the selected measures for this trial and time point
                for iC = 1:Ncomnds;
                    try
                        eval(thisfunCommands{iC});
                    catch
                        keyboard
                    end
                end
                
                clear temp thisPhi thisDphi;
                
            end
            
            %Interpolate through upsampling if required:
            if strcmpi(cfg.upsample,'yes')
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...upsample from Ncalc points to Ncut:
                        PC.(MeasNames{measInds(iM),iSubM}) = resample( PC.(MeasNames{measInds(iM),iSubM}), cfg.Ncut, cfg.Ncalc );
                    end
                end
                cfg.Nout(1) = cfg.Ncut;
                cfg.timeOut = cfg.timeCut;
            end
            
            %Smooth through convolution with a window function if required:
            if isfield(cfg,'smoothWin')
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...smooth:
                        PC.(MeasNames{measInds(iM),iSubM}) = conv( PC.(MeasNames{measInds(iM),iSubM}), cfg.smoothWin, 'same' );
                    end
                end
            end
            
            
            
    end
    
    DphiCell{nm(1),nm(2)}=Dphi;
    PCcell{nm(1),nm(2)}=PC;
    
    clear PC Dphi;
    
    
    toc
    
end


%For the statistics..

%Initialize statistics' results' structure
statsRes = cell(size(PCcell));
    
%In case surrogate statistics are to be made...
if ~isempty(cfg.stats)
    
    %...initialize phi...
    phi=phi0;
    clear phi0;
    
    %.unload parameters...
    pointStatMethod=cfg.stats.pointStatMethod;
    multiStatfun=cfg.stats.multiStatfun;
    alpha=cfg.stats.alpha;
    tail=cfg.stats.tail;
    corrMultComp=cfg.stats.corrMultComp;
    
    %...initialize a cell to store the measures structures for the
    %surrogates...
    PCsurr=cell(1,cfg.stats.Nperm);
    
    %...and a configurations' structure, identical to cfg, except for
    %the surrogate statistics...
    cfgSurr = cfg;
    cfgSurr.stats=[];
    
    %...for each pair of the surrogate time series...
    for iST=1:cfg.stats.Nperm;
        
        %...create it...
        phiSurr = cfg.stats.surrfun(phi,2,cfg.Ntr);
        
        %...and calculate the measures with exactly the same setting
        %except for the surrogate statistics...
        disp(['...calculating measures for surrogate data set ',num2str(iST),'/',num2str(cfg.stats.Nperm),'...'])
        [PCsurr{iST}, dummy, dummy2, dummy3] = DPphaseCoupling(phiSurr, cfgSurr, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas);
        
    end
    
    
    %...for each n,m combination...
    for iNM=1:K;
    
        %...this nm...:
        nm=cfg.nm(iNM,:);
        disp(['Calculating statistics for ',num2str(nm(1)),':',num2str(nm(2)),' coupling...'])
        
        %...initialize statistics' results' structure...
        statsRes{nm(1),nm(2)} = struct();
        
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
                        thisPCsurr =  zeros(cfg.Ntr,cfg.stats.Nperm);
                        for iST=1:cfg.stats.Nperm;    
                            %...for each permutation...
                            thisPCsurr(:,iST) = PCsurr{iST}{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM});
                        end     
                        
                        %...calculate the specific statistic for each trial...
                        [pointStat,  dummy] = DPcalcStat(PCcell{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}),thisPCsurr,pointStatMethod, multiStatfun);
                        clear thisPCsurr;
                        
                        %...calculate significance for each trial...
                        [statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).p,...
                         statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pUnCorr,...
                         statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).sign ] = ...
                         ...
                         DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                        
                        %...store statistics except for the ones of the surrogates...
                        %pointStat = rmfield(pointStat,'surrVal');
                        statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).stat = pointStat;
                        
                        
                        %...if trials are more than one do the same for the trials' mean...
                        if cfg.Ntr>1
                                                 
                            %...create a temporary matrix of the measure of the surrogate dataset...
                            thisPCmeanSurr =  zeros(1,cfg.stats.Nperm);
                            for iST=1:cfg.stats.Nperm;
                                thisPCmeanSurr(iST) = PCsurr{iST}{nm(1),nm(2)}.trialMean.(MeasNames{measInds(iM),iSubM});
                            end
                            
                            %...calculate the specific statistic...
                            [pointStat,  dummy] = DPcalcStat(PCcell{nm(1),nm(2)}.trialMean.(MeasNames{measInds(iM),iSubM}),thisPCmeanSurr,pointStatMethod, multiStatfun);
                            clear thisPCmeanSurr;
                            
                            %...and calculate significance...
                            [ statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMean,...
                              statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMeanUnCorr,...
                              statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).signMean ] = ...
                              ...
                              DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                          
                          %...store statistics except for the ones of the surrogates...
                          %pointStat = rmfield(pointStat,'surrVal');
                          statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).meanStat = pointStat;
                          
                        end
                        

                        
                    case 2 %trialTime
                        
                        %...create a temporary matrix of the measure of the surrogate dataset...
                        thisPCsurr =  zeros(cfg.Nout(1),cfg.Ntr,cfg.stats.Nperm);
                        for iST=1:cfg.stats.Nperm;
                            %...for each trial...
                            thisPCsurr(:,:,iST) = PCsurr{iST}{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM});
                        end
                        
                        %...calculate the specific statistic for each trial...
                        [pointStat,  multiStat] = DPcalcStat(PCcell{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}),thisPCsurr,pointStatMethod, multiStatfun);
                        clear thisPCsurr;
                        
                        %...calculate significance for each trial for each time point...
                        [statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).p,...
                         statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pUnCorr,...
                         statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).sign ] = ...
                         ...
                         DPcalcSignif(pointStat,alpha,tail,corrMultComp);

                        %...store statistics except for the ones of the surrogates...
                        %pointStat = rmfield(pointStat,'surrVal');
                        statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pointStat = pointStat;
                        
                        
                        if ~isempty(multiStat)
                            
                            %...calculate significance for each trial for multivariate statistic (all times)...
                            [statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMulti,...
                                statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMultiUnCorr,...
                                statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).signMulti ] = ...
                                ...
                                DPcalcSignif(multiStat,alpha,tail,corrMultComp);
                            
                            %...store statistics except for the ones of the surrogates...
                            %multiStat = rmfield(multiStat,'surrVal');
                            statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).multiStat = multiStat;
                        end
                        
                        %...if trials are more than one do the same for the trials' mean...
                        if cfg.Ntr>1
                            
                            %...create a temporary matrix of the measure of the surrogate dataset...
                            thisPCmeanSurr =  zeros(cfg.Nout(1),cfg.stats.Nperm);
                            for iST=1:cfg.stats.Nperm;
                                thisPCmeanSurr(:,iST) = PCsurr{iST}{nm(1),nm(2)}.trialMean.(MeasNames{measInds(iM),iSubM});
                            end
                            
                            %...calculate the specific statistic...
                            [pointStat,  multiStat] = DPcalcStat(PCcell{nm(1),nm(2)}.trialMean.(MeasNames{measInds(iM),iSubM}),thisPCmeanSurr,pointStatMethod, multiStatfun);
                            clear thisPCmeanSurr;
                            
                            %...and calculate significance for each time point....
                            [ statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMean,...
                              statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMeanUnCorr,...
                              statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).signMean ] = ...
                              ...
                              DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                          
                          %...store statistics except for the ones of the surrogates...
                          %pointStat = rmfield(pointStat,'surrVal');
                          statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).meanPointStat = pointStat;
                        
                          if ~isempty(multiStat)
                              
                              %...calculate significance for multivariate statistic (all times)...
                              [ statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMultiMean,...
                                  statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMultiMeanUnCorr,...
                                  statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).signMultiMean ] = ...
                                  ...
                                  DPcalcSignif(multiStat,alpha,tail,corrMultComp);
                              
                              %...store statistics except for the ones of the surrogates...
                              %multiStat = rmfield(multiStat,'surrVal');
                              statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).meanMultiStat = multiStat;
                              
                          end
                        end
                        
                        
                        
                    otherwise 
                        
                        %...create a temporary matrix of the measure of the surrogate dataset...
                        thisPCsurr =  zeros(cfg.Nout(1),cfg.stats.Nperm);
                        for iST=1:cfg.stats.Nperm;  
                            thisPCsurr(:,iST) = PCsurr{iST}{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM});
                        end
                        
                        %...calculate the specific statistic...
                        [pointStat,  multiStat] = ...
                         DPcalcStat(PCcell{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}),thisPCsurr,pointStatMethod, multiStatfun);
                        clear thisPCsurr;
                        
                        %...calculate significance for each time point...
                        [statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).p,...
                         statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pUnCorr,...
                         statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).sign ] = ...
                         ...
                         DPcalcSignif(pointStat,alpha,tail,corrMultComp);
                     
                        %...store statistics except for the ones of the surrogates...
                        %pointStat = rmfield(pointStat,'surrVal');
                        statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pointStat = pointStat;
                 
                        if ~isempty(multiStat)
                            %...calculate significance for multivariate statistic (all times)...
                            [statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMulti,...
                                statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).pMultiUnCorr,...
                                statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).signMulti ] = ...
                                ...
                                DPcalcSignif(multiStat,alpha,tail,corrMultComp);
                            
                            %...store statistics except for the ones of the surrogates...
                            %multiStat = rmfield(multiStat,'surrVal');
                            statsRes{nm(1),nm(2)}.(MeasNames{measInds(iM),iSubM}).multiStat = multiStat;
                        end
                         
                         
                end %method selection
                
                
                
                
            end %submeasure
        end %measure
        
        
        
    end %n,m combination
    
    
end


%Finish by extracting the results
Dphi=DphiCell;
PC=PCcell;

clear PCcell DphiCell;


disp('DONE!')



  


