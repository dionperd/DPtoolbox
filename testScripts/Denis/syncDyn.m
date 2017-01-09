%Toolbox locations: /Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox
%Put those into MATLAB path:
%io: input output between VA and MATLAB
%Time-Frequency: Gabor Transform
%Coupling Toolbox
%Statistics (if statistics are to be made)
%Utilities

clear all;
clc;

%A. Load VA exported Data into MATLAB
%(!export from VA as binary with format IEEE_FLOAT_32 and multiplexed!)
%Inputs:
filename = './DPdata.eeg';% full path, filename and extension of the VA data file
%Optional:
%channels: a vector of channel numbers [1 2 4 etc], or string 'all'
%trials: a vector of trials/segmetns numbers [1 2 4 etc], or string 'all'
[data hdr mrk stimulusMRK, responseMRK, segmentMRK, timezeroMRK] = DPreadBV(filename);
%data has and should have the format time x channels x trials/segments

%B. Gabor time-frequency transform
%Inputs:
%-x: data
%-f: the vector of frequencies 
%-Fs: sampling frequency in Hz
%-transform: either "gabor" or "gaborWavelet" (Morlet)
%Optional:
%-param: for Gabor transform is the gamma parameter, default: same as
%Gruber's software, i.e. gamma=2*sqrt(pi)
%[TF, param] = DPtfviaFFT(x,f,Fs,transform,param)
f=[4 6 8 10 12 20 30 40 60];%1 2 
[TF, dummy] = DPtfviaFFT(data,f,250,'Gabor'); 
%TF has the format time x frequency x channels x trials/segments and it is 
%complex transform

%C. Get the phase:
phiData=angle(TF);

%clear TF; %you don't need it anymore


%D. Parameters of configutation structure:

%1. Obligatory:
cfg.fs=hdr.Fs; %Sampling frequency in Hz
cfg.domain = 'phase'; %otherwise 'real' and a hilbert transform is performed
cfg.method = 'trialTime';%this means a sliding window within trial,
                         %other: 'trial' (1 value per trial, no time window), 'ensemble' (sliding window across trials)
%cfg.fc = [f(1) f(2)]; %central frequencies of signals, see below

%2. Optional but recommended:
cfg.measures = {'PLV', 'SE', 'CP','PLI','ICI'}; %default: 'all' (Notice 'IC' for all Viktor's measures together)
cfg.time = [-124:438]/cfg.fs;% in seconds, default: cfg.time = [0 N-1]/cfg.fs;
cfg.nm = [1 1]; %i.e., Dphi = n*phi1-m*phi2, several rows for testing different combinations
                %e.g., [1 1; 1 2; 2 1] etc, default = lcm(cfg.fc(1),cfg.fc(2))/cfg.fc
                %where lcm is the least common multiplier
%if domain=='real' and a hilbert transform is needed:
%   -cutTails: parts of the time series at the begining and end to be cut 
%                   out of the calculation in secs, vector of 2
%                   0=<positive real<N/fs/2 values, default: [TcMax TcMax]
%                   where TcMax=1/min(fc)

%3. In case that 'trialTime', or 'ensemble' is selected:
%The time points (sliding window centers) where the calculation will be
%performed. It has to be a subset of cfg.time. default: cfg.timeCalc =
%cfg.time, for pointwise calculation.
cfg.timeCalc = cfg.time; 


%4. More optional parameters:
%   -phaseTrnsfrm: 'yes' or 'no',string, default='no' 
%                but only if (Ncut or winlen)>TcMax (look below) 
% see Kralemann B, Cimponeriu L, Rosenblum MG, Pikovsky A, Mrowka R. 2008.
% Phase dynamics of coupled oscillators reconstructed from data. Phys Rev E. 77:1?16.
%
%   -Nbins: number of bins, positive real odd integer, 
%           default: calculated from data  
%           in the interval [5 33] as 2*floor(sqrt(N)/2)+1

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

%5. Parameters of measures:
%cfg.IC.Dphi0 = pi/4; (default)


%E. Optionally if statistics are required:
%cfg.stats.method = 'trialshuffling'; %Method of surrogates: 
%                                 Oher options: % -"timeshuffling": random time shuffling 
%                                                 -"phaseshuffling": random phase shuffling in Fourier space
%                                                 -"phaserandom": randomized phase in Fourier space
%                                                 default: "trialshuffling"
%cfg.stats.alpha = 0.05; %default = 0.01
%cfg.stats.tail =1; %default =1, otherwise, =2 for two tailed test
%Other options:
%       -Nperm: number of permutations, positive integer,
%               equal or greater than the default value: tail*ceil(1/alpha)+1
%       -pointStatMethod: a string defining the point statistic to be used,
%                         among 't' (for t-values),
%                               'z' for z-score normalization, or '' (none),
%                         default = '' 
%      Optionally, in the case of multiple comparisons correction
%      and if pointStat exists:
%       -corrMultComp: 'BONF' or 'FDR', default:'' for no correction
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


%F. calculation loop for all frequencies -not CFC, i.e., same frequency- and electrode pairs
cfg.Nch = 60; %60;
cfg.Nc = cfg.Nch*(cfg.Nch-1)/2;
cfg.Nf = length(f);%2;
cfg.f=f(1:cfg.Nf);

%For each measure we calculated...
for iM = 1:Nmeasures;
    %...and for each of the sumbmeasures of this measure ...
    for iSubM = 1:NmeasPmeas(measInds(iM));
        %...calculate the mean value (excluding possible nan
        %values)
        Sync.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,cfg.Ntr,cfg.Nc,cfg.Nf);
        %Sync.trialMean.(MeasNames{measInds(iM),iSubM})(:,:,iF,iC) = =zeros(cfg.Ncalc,cfg.Nc,cfg.Nf); 
    end
end


disp('Calculating synchronization...')
for iF=1:cfg.Nf; %for all frequencies
    
    cfg.winLen = 5/cfg.fs; %the window length in seconds, it will always have an odd number of points
    
    %Parameters
    cfg.fc = [cfg.f(iF) cfg.f(iF)]; %set the obligatory central frequencies vector
     
    disp(['fc = [',num2str(cfg.f(iF)),',',num2str(cfg.f(iF)),']'])
    
    iC=0;
    %for all channel pairs
    for iC1=1:cfg.Nch;
        
        for iC2=iC1+1:cfg.Nch; %No self couplings are included

            iC = iC+1;
            
            disp(['Channels = [',num2str(iC1),',',num2str(iC2),']'])
            
            %[phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds] = DPphaseCouplingPrepare(x,cfg);
            %Inputs:
            x(:,1,:) = squeeze(phiData(:,iF,iC1,:));
            x(:,2,:) = squeeze(phiData(:,iF,iC2,:));
            [phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds] = DPphaseCouplingPrepare(x,cfg);
            clear x;
            [PC, Dphi, cfg, THISstatsRes, MeasNames, NmeasPmeas] = DPphaseCoupling(phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds);
            
            %The result is given in a cell of structures in the form:
            %PC{n,m}.Measurename  = time (window) x trial matrix
            

            %For each measure we calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure ...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...calculate the mean value (excluding possible nan
                    %values)
                    Sync.(MeasNames{measInds(iM),iSubM})(:,:,iC,iF) = PC.(MeasNames{measInds(iM),iSubM});
                    %Sync.trialMean.(MeasNames{measInds(iM),iSubM})(:,:,iC,iF) = PC.trialMean.(MeasNames{measInds(iM),iSubM});
                    
                end
            end
            
          
            %Accordingly for the results of the statistics, one finds inside:
            %p values,
            %p_UnCorr for NO correction for multiple comparisons
            %sign: a matric of 1s (or 0s) for the values (time points x trials) that are
            %(not) significant
            % and some of the statistics of the surrogate data set
            
            
           % statsRes{iF,iC} = THISstatsRes{1,1};
            
            %......for the rest of the measures
        end
    end
end
save(['syncDyn',cfg.method,'.mat'],'PLV','CP','SE','cfg');%,'PCI','NCI','ACI','ICI12','ICI21','mPLV','mPCI','mNCI','mACI','mICI12','mICI21','cfg','statsRes')

load('syncDyntrialTime.mat')


%For each measure we calculated...
for iM = 1:Nmeasures;
    %...and for each of the sumbmeasures of this measure ...
    for iSubM = 1:NmeasPmeas(measInds(iM));
        %...calculate the mean value (excluding possible nan
        %values)
        Corr.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,cfg.Ncalc,cfg.Ntr,cfg.Nf);
        PCA.coefs.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Nc,cfg.Nc,cfg.Ntr,cfg.Nf);
        PCA.scores.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Ncalc,cfg.Nc,cfg.Ntr,cfg.Nf);
        PCA.vars.(MeasNames{measInds(iM),iSubM})=zeros(cfg.Nc,cfg.Ntr,cfg.Nf);
    end
end



for iF=1:cfg.Nf;
   
    for iTr=1:cfg.Ntr;
        
        disp('')
        disp(['f = ',num2str(cfg.f(iF))])
        disp(['trial = ',num2str(iTr)])
        
        tic
        disp('Calculating correlations...')
        for iT1 = 125:100:325;
            
            for iT2 = 125:100:425
                
                disp(['Times = [',num2str(iT1),',',num2str(iT2),']'])
                
                %For each measure we calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure ...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...calculate the mean value (excluding possible nan
                        %values)
                        C1 = squeeze( Sync.(MeasNames{measInds(iM),iSubM})(iT1,iTr,:,iF) );
                        C2 = squeeze(Sync.(MeasNames{measInds(iM),iSubM})(iT2,iTr,:,iF) );
                        Corr.(MeasNames{measInds(iM),iSubM}) = xcorr(C1,C2,0,'coeff');
                    end
                end
                
            end
        end
        toc
        
        disp('...and pca...')
        tic
        %For each measure we calculated...
        for iM = 1:Nmeasures;
            %...and for each of the sumbmeasures of this measure ...
            for iSubM = 1:NmeasPmeas(measInds(iM));
                %...calculate the mean value (excluding possible nan
                %values)
                tempC = squeeze( Sync.(MeasNames{measInds(iM),iSubM})(:,iTr,:,iF) );
                [PCA.coefs.(MeasNames{measInds(iM),iSubM})(:,:,iTr,iF),...
                 PCA.scores.(MeasNames{measInds(iM),iSubM})(:,:,iTr,iF),...
                 PCA.vars.(MeasNames{measInds(iM),iSubM})(:,iTr,iF)]...
                                                  = princomp(zscore(tempC));
                PCA.vars.(MeasNames{measInds(iM),iSubM})(:,iTr,iF) = ...
                 PCA.vars.(MeasNames{measInds(iM),iSubM})(:,iTr,iF) / ...
                 sum(PCA.vars.(MeasNames{measInds(iM),iSubM})(:,iTr,iF)); %calculate percentage of variance explained
            end
        end

        toc
        
    end
    
end


save(['syncDyn',cfg.method,'.mat'],'Corr','PCA','-append');