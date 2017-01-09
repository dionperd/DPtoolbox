function phaseCoupling_start(filenum)

cd('/Volumes/LIP-Git2/Caroline/Gaming_study_downgesampled/Export_down');
if ischar(filenum)
filenum=str2num(filenum);
end

FolderDir=dir;
files=cell(2,1);
ii=0;
for iFile = 3:length(dir);
    if any(strcmpi(FolderDir(iFile).name(end-2:end),{'dat'}));%{'eeg','dat','bin'}
        ii=ii+1;
        files{ii}=FolderDir(iFile).name;
    end
end

%%% ...and run the CM-calculation for Input-specified files
for iFile = filenum;
    disp('-----------------')
    disp(['TORTURE FILE ',files{iFile}])
    
    filename=files{iFile}


%Toolbox locations: /Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox
%Also in Dropbox now!
%Put those into MATLAB path:
%io: input output between VA and MATLAB
%Time-Frequency: Gabor Transform
%Coupling Toolbox
%Statistics (if statistics are to be made)

% clear all;
% clc;

%A. Load VA exported Data into MATLAB
%(!export from VA as binary with format IEEE_FLOAT_32 and multiplexed!)
%Inputs:
%filename = '/Volumes/LIP-Git2/Caroline/Gaming_study_downgesampled/Export_down.eeg';% full path, filename and extension of the VA data file
%Optional:
%channels: a vector of channel numbers [1 2 4 etc], or string 'all'

channels=[1;2;3;8;10;12;14;16;26;28;30;32;34;44;46;48;50;52;58;59;60;65;66;67;72;74;76;78;80;90;92;96;94;98;108;110;112;114;116;122;123;124];
%trials: a vector of trials/segmetns numbers [1 2 4 etc], or string 'all'
[data, hdr, mrk, dummy] = DPreadBV(filename,channels);
%data has and should have the format time x channels x trials/segments
[N, Nch, Ntr] = size(data);
cfg.Nch = Nch;

%For the 2 Hz bin calculation:
%The centers of the frequency bands
f=[2 4 6 8 10 12 14 16 18 20];
cfg.Nf=length(f);


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
%[TF, dummy] = DPtfviaFFT(data,f,250,'Gabor'); 
c=5; %wavelet parameter: number of cycles per time window for each frequency
[TF, dummy] = DPtfviaFFT(data,f,hdr.Fs,'gaborWavelet',c);
%TF has the format time x frequency x channels x trials/segments and it is 
%complex transform


%C. Get the phase:
phiData=angle(TF);
clear TF; %you don't need it anymore


%D. Parameters of configutation structure:

%1. Obligatory:
cfg.fs=hdr.Fs; %Sampling frequency in Hz
%For the broad frequency bands calculation:
%cfg.domain = 'real';
cfg.domain = 'phase'; %otherwise 'real' and a hilbert transform is performed
cfg.method = 'trialTime';%this means a sliding window within trial,
                         %other: 'trial' (1 value per trial, no time window), 'ensemble' (sliding window across trials)
%cfg.fc = [f(1) f(2)]; %central frequencies of signals, see below

%2. Optional but recommended:
cfg.measures = {'PLV', 'CCR','IC'}; %default: 'all' (Notice 'IC' for all Viktor's measures together)
cfg.time = [1:size(phiData,1)]/cfg.fs;% in seconds, default: cfg.time = [0 N-1]/cfg.fs;
cfg.nm = [1 1]; %i.e., Dphi = n*phi1-m*phi2, several rows for testing different combinations
                %e.g., [1 1; 1 2; 2 1] etc, default = lcm(cfg.fc(1),cfg.fc(2))/cfg.fc
                %where lcm is the least common multiplier
%if domain=='real' and a hilbert transform is needed:
cfg.cutTails=[0 0]; %For the broad frequency bands calculation:
%   -cutTails: parts of the time series at the begining and end to be cut 
%                   out of the calculation in secs, vector of 2
%                   0=<positive real<N/fs/2 values, default: [TcMax TcMax]
%                   where TcMax=1/min(fc)

%3. In case that 'trialTime', or 'ensemble' is selected:
%The time points (sliding window centers) where the calculation will be
%performed. It has to be a subset of cfg.time. default: cfg.timeCalc =
%cfg.time, for pointwise calculation.
cfg.timeCalc = cfg.time(1:25:end); %one point every 25, i.e., one calculation at each 100 ms, given that Ts = 4ms or fs=250Hz.
Ncalc = length(cfg.timeCalc);
cfg.winLen = 0.8; %the window length in seconds, it will always have an odd number of points

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
%when timeCalc is a downsampled version of time or timeCut, then
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
cfg.IC.Dphi0 = pi/4; %(default)


%E. Optionally if statistics are required:
cfg.stats=[];
%cfg.stats.method = 'trialshuffling'; %Method of surrogates: 
%                                 Oher options: % -"timeshuffling": random time shuffling 
%                                                 -"phaseshuffling": random phase shuffling in Fourier space
%                                                 -"phaserandom": randomized phase in Fourier space
%                                                 default: "trialshuffling"
%cfg.stats.alpha = 0.05; %default = 0.05
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
OutMeasures = {'PLV','CCR','ICI12','ICI21'};
NoutMeasures = length(OutMeasures);

tic
%Initialize loop
Nloops = cfg.Nf*(cfg.Nch*(cfg.Nch-1)/2);
INDs = nan(Nloops,3);
fc = nan(Nloops,1);
Labels = cell(Nloops,3);
iL=0;
for iF=1:cfg.Nf; %for all frequencies
    
    for iC1=1:cfg.Nch; %for all channel pairs
        
        for iC2=iC1+1:cfg.Nch; %No self couplings are allowed
            
            iL=iL+1;
            
            INDs(iL,:) = [iF iC1 iC2]; %the indices of center frequencies and channels
            
            fc(iL) = f(iF); %the central frequency
            
            Labels(iL,:) = {num2str(f(iF)),hdr.label{iC1},hdr.label{iC2}}; %frequency and channel labels
        end
    end
end
toc

%For each measure we will calculate...
for iM = 1:NoutMeasures;
    %...calculate the mean value (excluding possible nan
    %values)
    C.(OutMeasures{iM})=zeros(Ncalc,Ntr,Nloops);
    %Cmean.(OutMeasures{iM})=zeros(Ncalc,Nloops);  
end


%Loop over all of the above calculations:
for iL=1:Nloops; 
    
    %Parameters
    cfg.fc = [fc(iL) fc(iL)]; %set the obligatory central frequencies vector    
    
    %[phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas] = DPphaseCouplingPrepare(x,cfg);
    %Inputs:
    x(:,1,:) = squeeze(phiData(:,INDs(iL,1),INDs(iL,2),:));
    x(:,2,:) = squeeze(phiData(:,INDs(iL,1),INDs(iL,3),:));
    [phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas] = DPphaseCouplingPrepare(x,cfg);
    clear x;
    [PC, Dphi, cfg, THISstatsRes] = DPphaseCoupling(phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds, MeasNames, NmeasPmeas);
    
    
    %The result is given in a cell of structures in the form:
    %PC{n,m}.Measurename  = time (window) x trial matrix
    
    for iM = 1:NoutMeasures;
        %...calculate the mean value (excluding possible nan
        %values)
        C.(OutMeasures{iM})(:,:,iL)=PC{cfg.nm(1),cfg.nm(2)}.(OutMeasures{iM});
        %Cmean.(OutMeasures{iM})(:,iL) = PC{cfg.nm(1),cfg.nm(2)}.trialMean.(OutMeasures{iM});
    end

%     %For each measure we calculated...
%     for iM = 1:Nmeasures;
%         %...and for each of the sumbmeasures of this measure ...
%         for iSubM = 1:NmeasPmeas(measInds(iM));
%             %...calculate the mean value (excluding possible nan
%             %values)
%             C.(MeasNames{measInds(iM),iSubM})(:,:,iL)=PC{cfg.nm(1),cfg.nm(2)}.(MeasNames{measInds(iM),iSubM});
%             %Cmean.(OutMeasures{iM})(:,iL) = PC{cfg.nm(1),cfg.nm(2)}.trialMean.(MeasNames{measInds(iM),iSubM});
%         end
%     end
    
    %......for the rest of the measures
    
end
cfg.fc=fc;
cfg.nm=nm;
cfg.Labels = Labels;
[pathstr, filename, ext] = fileparts(filename);
save([filename,'_PC.mat'],'C','cfg');%'Cmean',

end
