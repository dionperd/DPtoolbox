function DPdata = DPcreateData(domain,fsample,time,data,varargin)%(tau/freq),datasetLabels,fileName,filePath,saveFlag,axesLabels,axesUnits,axesScales,axesValues

%This function creates a DPdata structure for time series data.

%DPdata is the main data structure creator of the DP toolbox

%signals: 2D cell array of tstool 2D signal objects that contains the data 
%      Form and indexing:
%      signals{i,j} is indexing a single dataset (signal (e.g. channel or component) x trial):
%       
%      data{..} indexing goes as follows:

%      Optionally only for time frequency data:
%      signals is a cell of tstool 3D signal objects with size(dataset,3)=2
%                
%      each signals{...} element should be arranged columnwise,
%      whereas it may be empty, as well
%
% datasetLabels: cell array of strings of numel equal to numel(size(signals)). 
%        Each string corresponds to the label of each index of data{...}
%        Default is {'',''}

% domain: string determining the domain of the data ('time','time-freq', or
% 'state space')

% fsample: scalar for sampling frequency (in Hz) assumed to be the same 
%          for all single datasets
%

% time: a column vector of monotonously increasing values that accounts for the
%       time vector to which the data refers

% dims:  array (dims(i,j)) of the dimensions of datasets, dims=0 for
%        empty datasets, 1 for time series data, length(freqs) for time-frequency data
%        and >=1 for state space data%   

% tau: cell array (tau{i,j}) of vectors of nonnegative values
%              determining the time delay of the mth time-delay vector 
%              from the (m-1)th one  

% freq: a row vector of monotonously increasing values that accounts for the
%       frequency vector to which the data refers

% axes: a cell array containing tstool achses objects, as many as the axes
%       of the data (time-signal for time data, datapoint-signal for state
%       space data and time-frequency-real part/amplitude/squared amplitude
%       and imaginative part/phase ofr time=frequency data)

%
% file: a structure containing the following fields
%      -name: string for the file name of the dataset WITH the proper extension
%      -path: string for file path of the dataset 
%
%      This is not an input argument-the creator function calculates it:
%      -size: double for the size of the file in MBs
% %
% 
% 
% INPUTS:
% 
% Obligatory:
% -fsample
% -time
% -domain
% -data 
% -tau/freq only for state space/time-frequency data respectively

% Optional:
% -datasetLabels
% -fileName: a string of the file name (with extension) where the DPdata should be saved
%            if it is not given, DPdata will not be saved
% -filePath: a string of the file path. if is is not given, DPdata will be
%            saved in the current path
% -saveFlag: "Yes" or "No" to save the structure or not
% -axesLabels: cell of 2 strings, each of which labelling an axis
% -axesUnits: cell of 2 strings, each of which determines the units of an axis
% -axesScales: cell of 2 strings, each of which determines the scale of an axis (it can be 'linear', 'log2','log10','ln' or 'arbitrary')
%              for linear (logarithmic) scaling, we assume that the
%              corresponding data are already linear (logarithmic)
% -axesValues: cell of 2 row vectors of doubles, determining the values of each axis that is signaled as of 'arbitrary'scale


% 
% OUTPUT:
% -DPdata structure



%Validate the inputs
%[RESULT x] = DPvalidateData(x,testfun,param,mode,execfun,default,varName,funcName)
funcName='DPcreateData';

%Required inputs
varName='domain';
testDOMAIN = { @(domain)ischar(domain),...
               @(domain)isvector(domain),... 
               @(domain)any( strcmpi(domain,'time')||strcmpi(domain,'time-freq')||strcmpi(domain,'state space')) };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default=nan;
[RESULT domain] = DPvalidateData(domain,testDOMAIN,param,mode,execfun,default,varName,funcName);
%Put domain inside structure
DPdata.domain=domain;

varName='fsample';
testFSAMPLE = {@(fsample)isnumeric(fsample),...
               @(fsample)isscalar(fsample),...
               @(fsample)(~isempty(fsample)),...
               @(fsample)isreal(fsample),...
               @(fsample)( fsample>0 ),...
               @(fsample)isa(fsample,'double') };
param={{},{},{},{},{},{}};
mode=['e','e','e','e','e','w'];
execfun={{},{},{},{},{},@(time)double(time)};
default=nan;
[RESULT fsample] = DPvalidateData(fsample,testFSAMPLE,param,mode,execfun,default,varName,funcName);
%Put fsample into DPdata structure
DPdata.fsample=fsample;

varName='time';
testTIME = { @(time)isnumeric(time),...
             @(time)isvector(time),... 
             @(time)isreal(time),...
             @(time) all(diff(time)>0)...,
             @(time)(size(time,1)>=size(time,2)),...
             @(time)isa(time,'double')};
param={{},{},{},{},{},{}};
mode=['e','e','e','e','w','w'];
execfun={{},{},{},{},@(time)time.',@(time)double(time)};
default=nan;
[RESULT time] = DPvalidateData(time,testTIME,param,mode,execfun,default,varName,funcName);
Ntimes=length(time);
%Put time into DPdata structure
DPdata.time=time;

varName='data';
testDATA = { @(data)iscell(data),...
             @(data)(~isempty(data)),...
             @(data)length(size(data))==2};
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default=nan;
[RESULT data] = DPvalidateData(data,testDATA,param,mode,execfun,default,varName,funcName);

N_datasets = numel(data); %the number of datasets
dataSIZE=size(data);%the size of data cell
    
%Initialize dims
dims=zeros(dataSIZE);

if strcmpi(domain,'time')
    
    OPT_ARG = 0;
    
    Naxes=2;
    defaultAxesLabels={'time',''};
    defaultAxesUnits={'s',''};
    
    
    for ii=1:N_datasets;
        
        varName=['data{',num2str(ii),'}'];
        
        testDATASET = {@(dataset)isnumeric(dataset),...
                       @(dataset)( length(size(dataset))==2 ),...
                       @(dataset)isreal(dataset),...
                       @(dataset)isa(dataset,'double'),...
                       @(dataset,Ntimes) isempty(dataset)||any( Ntimes==size(dataset) ),...
                       @(dataset,Ntimes) isempty(dataset)||( size(dataset,1)==Ntimes ) };
        
        param={{},{},{},{},{Ntimes},{Ntimes}};
        mode=['e','e','e','w','e','w'];
        execfun={{},{},{},@(dataset)double(dataset),{}, @(dataset)dataset.'};
        default=nan;
        [RESULT data{ii}] = DPvalidateData(data{ii},testDATASET,param,mode,execfun,default,varName,funcName);
        
        %...store dimensions of the dataset...
        dims(ii)=size(data{ii},2);
        
    end
    

elseif strcmpi(domain,'time-freq')
    
    if dataSIZE(end)~=2
        error('dataSIZE(end)~=2 for time-frequency data')
    end
      
    OPT_ARG = 1;
    
    Naxes=4;
    defaultAxesLabels={'time','frequency','real','imag'};
    defaultAxesUnits={'s','Hz','',''};
    
    freq=varargin{1};
    
    varName='freq';
    testFREQ = { @(freq)isnumeric(freq),...
                 @(freq)isvector(freq),...
                 @(freq)isreal(freq),...
                 @(freq)all(freq>0)...,
                 @(freq)all(diff(freq)>0)...,
                 @(freq)(size(freq,2)>=size(freq,1)),...
                 @(freq)isa(freq,'double')};
    param={{},{},{},{},{},{},{}};
    mode=['e','e','e','e','e','w','w'];
    execfun={{},{},{},{},{},@(freq)freq.',@(freq)double(freq)};
    default=nan;
    [RESULT freq] = DPvalidateData(freq,testFREQ,param,mode,execfun,default,varName,funcName);
    Nfreqs=length(freq);
    %Put freq inside the DPdata structure
    DPdata.freq=freq;
    DPdata.dims=Nfreqs;
    
    for ii=1:N_datasets;
        
        varName=['data{',num2str(ii),'}'];
        
        
        datasetSIZE = size(data{ii});
        
        testDATASET = {@(dataset)isnumeric(dataset),...
                       @(dataset)( length(size(dataset))==3 ),...
                       @(dataset)( ~isreal(dataset) ),...
                       @(dataset,Ntimes,Nfreqs,datasetSIZE)(...
                         ( datasetSIZE(1)==Ntimes && datasetSIZE(2)==Nfreqs && datasetSIZE(3)==2 ) ||...
                         ( datasetSIZE(1)==Nfreqs && datasetSIZE(2)==Ntimes && datasetSIZE(3)==2) ||...
                          all(datasetSIZE==0)                                 ),...
                       @(dataset,SIZES)(...
                         all(datasetSIZE==0) ||...
                        ( datasetSIZE(1)==Ntimes && datasetSIZE(2)==Nfreqs ) ),...
                       @(dataset)isa(dataset,'double')};
        
        
        param={{},{},{},{Ntimes,Nfreqs,datasetSIZE},{Ntimes,Nfreqs,datasetSIZE},{}};
        mode=['e','e','e','e','w','w'];
        execfun={{},{},{},{}, @(dataset)permute(dataset,[2,1,3]),@(dataset)double(dataset),};
        default=nan;
        [RESULT data{ii}] = DPvalidateData(data{ii},testDATASET,param,mode,execfun,default,varName,funcName);
        

    end
        
else
    
    OPT_ARG = 1;
    
    Naxes=2;
    defaultAxesLabels={'point',''};
    defaultAxesUnits={'',''};
    
    tau=varargin{1};
    
    varName='tau';
    testTAU = { @(tau)iscell(tau),...
                @(tau,dataSIZE)size(tau)==dataSIZE };
    param={{},{dataSIZE}};
    mode=['e','e'];
    execfun={{},{}};
    default=nan;
    [RESULT tau] = DPvalidateData(tau,testTAU,param,mode,execfun,default,varName,funcName);
    

    %validate each separate dataset and the respective tau
    for ii=1:N_datasets;
        
        varName=['tau{',num2str(ii),'}'];
        
        testTAU_i = { @(tau_i)isnumeric(tau_i),...
                      @(tau_i)( length(size(tau_i))==2 ),...
                      @(tau_i)isreal(tau_i),...
                      @(tau_i)all(tau_i>=0),...
                      @(tau_i)all(round(tau_i)==tau_i),...
                      @(tau_i)( size(tau_i,1)<=size(tau_i,2) ),...
                      @(tau_i)isa(tau_i,'double') };
        
        param={{},{},{},{},{},{},{}};
        mode=['e','e','e','e','e','w','w'];
        execfun={{},{},{},{},{},@(tau_i)tau_i.', @(tau_i)double(tau_i)};
        default=nan;
        [RESULT tau{ii}] = DPvalidateData(tau{ii},testTAU_i,param,mode,execfun,default,varName,funcName);
        
        tauWIN = sum(tau{ii});
        
        Npoints=Ntimes-tauWIN;
        
        
        %...store dimensions of the dataset...
        dims(ii)=length(size(tau{ii}));
        
        
        varName=['data{',num2str(ii),'}'];
        
        datasetSIZE = size(data{ii});
        
        testDATASET = {@(dataset)isnumeric(dataset),...
                       @(dataset)( length(size(dataset))==2 ),...
                       @(dataset)isreal(dataset),...
                       @(dataset,Ntimes,Npoints,dim,datasetSIZE)(...
                        ( any(datasetSIZE(1)==[Ntimes,Npoints]) && datasetSIZE(2)==dim ) ||...
                        ( datasetSIZE(1)==dim && any(datasetSIZE(2)==[Ntimes,Npoints]) ) ||...
                         all(datasetSIZE==0) && (dim==0) ),...
                       @(dataset,Ntimes,N,dim,datasetSIZE)(...
                         all(datasetSIZE==0) && (dim==0) ||...
                        ( any(datasetSIZE(1)==[Ntimes,Npoints]) && datasetSIZE(2)==dim ) ),...
                       @(dataset)isa(dataset,'double') };
        
        param={{},{},{},{Ntimes,Npoints,dims(ii),datasetSIZE},{Ntimes,Npoints,dims(ii),datasetSIZE},{}};
        mode=['e','e','e','e','w','w'];
        execfun={{},{},{},{}, @(dataset,Ntimes,N,dim,datasetSIZE)dataset.', @(dataset)double(dataset)};
        default=nan;
        [RESULT data{ii}] = DPvalidateData(data{ii},testDATASET,param,mode,execfun,default,varName,funcName);
        
        %Add nans so that all datasets have the same number of data points
        %(same 1st dimension)
        if (tauWIN>0)&&(datasetSIZE(1)==Npoints)
            data{ii} = [data{ii};nan(tauWIN,dims(ii))];
        end
        
    end
   
       
    %Put tau inside the DPdata structure
    DPdata.tau=tau;
    
end
%Put dims inside the DPdata structure
DPdata.dims=dims;


%Unload optional inputs of varargin...
datasetLabelsDEFAULT={'',''};
if nargin>4+OPT_ARG
    datasetLabels=varargin{1+OPT_ARG};
        
    varName='datasetLabels';
    testLABELS = { @(datasetLabels)iscellstr(datasetLabels),...
                   @(datasetLabels)( numel(datasetLabels)==2 ) };
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    default=datasetLabelsDEFAULT;
    [RESULT datasetLabels] = DPvalidateData(datasetLabels,testLABELS,param,mode,execfun,default,varName,funcName);
    
    if RESULT==1
        
        for ii=1:2;
            
            varName=['datasetLabels{',num2str(ii),'}'];
            testLABELSn = { @(label)ischar(label),...
                            @(label)isvector(label)||strcmpi(label,'') };
            param={{},{}};
            mode=['e','e'];
            execfun={{},{}};
            default='';
            [RESULT datasetLabels{1}] = DPvalidateData(datasetLabels{1},testLABELSn,param,mode,execfun,default,varName,funcName);
            
        end
        
    end
    
else
    datasetLabels=datasetLabelsDEFAULT;
end

%Put datasetLabels into DPdata structure
DPdata.datasetLabels=datasetLabels;

if nargin>5+OPT_ARG
    fileName=varargin{2+OPT_ARG};
    
    varName='fileName';
    testFILENAME = {@(fileName)ischar(fileName),...
                    @(fileName)isvector(fileName)||strcmpi(fileName,'') };
    param={{},{}};    
    mode=['e','e'];
    execfun={{},{}};
    default='';
    [RESULT fileName] = DPvalidateData(fileName,testFILENAME,param,mode,execfun,default,varName,funcName);

else
    fileName='';
end


if nargin>6+OPT_ARG
    filePath=varargin{3+OPT_ARG};
    
    varName='filePath';
    testFILEPATH = {@(filePath)ischar(filePath),...
                    @(filePath)isvector(filePath)||strcmpi(filePath,'') }; 
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    default='';
    [RESULT filePath] = DPvalidateData(filePath,testFILEPATH,param,mode,execfun,default,varName,funcName);

else
    filePath='';
end

if nargin>7+OPT_ARG
    saveFlag=varargin{4+OPT_ARG};
    
    varName='saveFlag';
    testSAVEFLAG = {@(saveFlag)ischar(saveFlag),...
                    @(saveFlag)isvector(saveFlag),...
                    @(saveFlag) any(strcmpi(saveFlag,{'Yes','No'}))}; 
    param={{},{},{}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default='No';
    [RESULT saveFlag] = DPvalidateData(saveFlag,testSAVEFLAG,param,mode,execfun,default,varName,funcName);

else
    saveFlag='No';
end


if nargin>8+OPT_ARG
    
    axesLabels=varargin{5+OPT_ARG};
    
    varName='axesLabels';
    
    testAXESLABELS = { @(axesLabels)iscellstr(axesLabels),...
                       @(axesLabels,Naxes)( numel(axesLabels)==Naxes ) };
    param={{},{Naxes}};
    mode=['e','e'];
    execfun={{},{}};
    [RESULT axesLabels] = DPvalidateData(axesLabels,testAXESLABELS,param,mode,execfun,defaultAxesLabels,varName,funcName);
    
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    for ii=1:Naxes;
        if RESULT==1
            varName=['axesLabels{',num2str(ii),'}'];
            testAXESLABEL_i = { @(axesLabel_i)ischar(axesLabel_i),...
                                @(axesLabel_i)isvector(axesLabel_i)||strcmpi(axesLabel_i,'') };
            
             [RESULT axesLabels{ii}] = DPvalidateData(axesLabels{ii},testAXESLABEL_i,param,mode,execfun,defaultAxesLabels{ii},varName,funcName);
            
        end
    end
    
    
else
    axesLabels=defaultAxesLabels;
end


if nargin>9+OPT_ARG
    
    axesUnits=varargin{6+OPT_ARG};
    
    varName='axesUnits';
    
    testAXESUNITS = { @(axesUnits)iscellstr(axesUnits),...
                      @(axesUnits,Naxes)( numel(axesUnits)==Naxes ) };
    param={{},{Naxes}};
    mode=['e','e'];
    execfun={{},{}};
    [RESULT axesUnits] = DPvalidateData(axesUnits,testAXESUNITS,param,mode,execfun,defaultAxesUnits,varName,funcName);
    
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    for ii=1:Naxes;
        if RESULT==1
            varName=['axesUnits{',num2str(ii),'}'];
            testAXESUNIT_i = { @(axesUnit_i)ischar(axesUnit_i),...
                               @(axesUnit_i)isvector(axesUnit_i)||strcmpi(axesUnit_i,'') };
            param={{},{}};
            mode=['e','e'];
            execfun={{},{}};
            [RESULT axesUnits{ii}] = DPvalidateData(axesUnits{ii},testAXESUNIT_i,param,mode,execfun,defaultAxesUnits{ii},varName,funcName);
            
        end
    end
    
    
else
    axesUnits=defaultAxesUnits;
end


defaultAxesScales={'linear'};
defaultAxesScales=repmat(defaultAxesScales,1,Naxes);
if nargin>10+OPT_ARG
    axesScales=varargin{7+OPT_ARG};
    
    varName='axesScales';
    testAXESSCALES = { @(axesScales)iscellstr(axesScales)
                       @(axesScales,Naxes)( numel(axesScales)==Naxes ) };
                   
    param={{},{Naxes}};
    mode=['e','e'];
    execfun={{},{}};
    [RESULT axesScales] = DPvalidateData(axesScales,testAXESSCALES,param,mode,execfun,defaultAxesScales,varName,funcName);
    
    if RESULT==1
        
        for ii=1:Naxes;
            
            varName=['axesScales{',num2str(ii),'}'];
            testAXESSCALE_i = { @(axesScale_i)strcmpi(axesScale_i,'linear')||strcmpi(axesScale_i,'log2')||strcmpi(axesScale_i,'log10')||strcmpi(axesScale_i,'ln')||strcmpi(axesScale_i,'arbitrary') };
            param={{}};
            mode=['e'];
            execfun={{}};
            [RESULT axesScales{ii}] = DPvalidateData(axesScales{ii},testAXESSCALE_i,param,mode,execfun,defaultAxesScales{ii},varName,funcName);
            
        end
        
    end
else
    axesScales=defaultAxesScales;
end

defaultAxesValues={[]};
defaultAxesValues = repmat(defaultAxesValues,1,Naxes);
if nargin>11+OPT_ARG
    
    axesValues=varargin{8+OPT_ARG};
    
    varName='axesValues';
    testAXESVALUES = { @(axesValues)iscell(axesValues),...
                       @(axesValues,Naxes)( numel(axesValues)==Naxes ) };
    param={{},{Naxes}};
    mode=['e','e'];
    execfun={{},{}};
    [RESULT axesValues] = DPvalidateData(axesValues,testAXESVALUES,param,mode,execfun,defaultAxesValues,varName,funcName);
    
    if RESULT==1
        
        for ii=1:Naxes;
            
            if strcmpi(axesScales{ii},'arbitrary')
                varName=['axesValues{',num2str(ii),'}'];
                testAXESVALUE_i = { @(axesValue_i)isnumeric(axesValue_i),...
                                    @(axesValue_i)isvector(axesValue_i)||isempty(axesValue_i),...
                                    @(axesValue_i)isreal(axesValue_i),...
                                    @(axesValue_i)all( diff(axesValue_i)>0 ),...
                                    @(axesValue_i)isa(axesValue_i,'double') };
                param={{},{},{},{},{}};
                mode=['e','e','e','e','w'];
                execfun={{},{},{},{},@(axesValues)double(axesValues)};
                [RESULT axesValues{ii}] = DPvalidateData(axesValues{ii},testAXESVALUE_i,param,mode,execfun,defaultAxesValues{ii},varName,funcName);

            else
                axesValues{ii}=defaultAxesValues{ii};
            end
        end
        
    end
    
else
    axesValues=defaultAxesValues;
    for ii=1:Naxes;
        if strcmpi(axesScales{ii},'arbitrary') && isempty(axesValues{ii})
            cprintf('TEXT',['\nWARNING in function ',funcName,':','\naxesValues{',num2str(ii),'} is empty although axesScales{',num2str(ii),'}=''arbitrary''.\nSetting default axesScales{',num2str(ii),'}=''linear''.\n'])
            axesScales{ii} = 'linear';
        end
       
    end
end




%Create time axis:

%...create unit...
timeAXIS=defineTimeFreqAxes(time,axesLabels(1),axesUnits(1),axesScales(1),axesValues(1));


%LOG_AXES = strcmpi(axesScales{2},'log2')||strcmpi(axesScales{2},'log10')||strcmpi(axesScales{2},'ln');
if strcmpi(domain,'time-freq')
    
    %Create frequency axis:
    
    freqAXIS=defineTimeFreqAxes(freq,axesLabels(2),axesUnits(2),axesScales(2),axesValues(2));

    
    
    [dataAXES dataMINs dataMAXs dataDELTAs]=initializeDataAxes(axesLabels(3:4),axesUnits(3:4),axesScales(3:4),axesValues(3:4));
    
  
    
    
    for ii=1:N_datasets;
        
        [currAXES dataMINs dataMAXs dataDELTAs]=defineCurrDataAxes(dataAXES,dataMINs, dataMAXs, dataDELTAs,{data{ii}(:,:,1), data{ii}(:,:,2)},axesScales(3:4));
        
        %...and create signal
        data{ii} = signal(data{ii},timeAXIS,freqAXIS,currAXES{1},currAXES{2});
        
    end
    
    dataAXES=updateDataAxes(2,dataAXES,dataMINs, dataMAXs, dataDELTAs,axesScales(3:4));
    
    DPdata.axes={timeAXIS, freqAXIS, dataAXES{1}, dataAXES{2}};
    
else
    for ii=1:N_datasets;
        
        [dataAXIS dataMIN dataMAX dataDELTA]=initializeDataAxes(axesLabels(2),axesUnits(2),axesScales(2),axesValues(2));
        
        [currAXIS dataMIN dataMAX dataDELTA]=defineCurrDataAxes(dataAXIS,dataMIN, dataMAX, dataDELTA, {data{ii}},axesScales(2));
        
        %...and create signal
        data{ii} = signal(data{ii},timeAXIS,currAXIS);
        
    end
    
    dataAXIS=updateDataAxes(dataAXIS,dataMIN, dataMAX, dataDELTA,axesScales(2));
    
    DPdata.axes={timeAXIS, dataAXIS{:}};
    
end

%Creat DPdata structure and put data inside
DPdata.signals=data;
clear data;


%Create file structure

if length(fileName)>4
    if strcmpi(domain,'time')
        DPdata.file.name=[fileName(1:end-4),'_tm',fileName(end-3:end)];
    elseif strcmpi(domain,'time-freq')
        DPdata.file.name=[fileName(1:end-4),'_tf',fileName(end-3:end)];
    else 
        DPdata.file.name=[fileName(1:end-4),'_st',fileName(end-3:end)];
    end
else
    DPdata.file.name='';
end
%....check if the file path exists...
if exist(filePath,'dir')~=7
    %...if not set the default one...
    DPdata.file.path='';
else
    DPdata.file.path=filePath;
end

%...set the default size
DPdata.file.size=0;
%...get the actual size...
dummy=whos('DPdata');
%...and set it
DPdata.file.size=dummy.bytes/1048576; %MBs=bytes/1024^2 = bytes/1048576


%Saving DPdata in the hard disk
if strcmpi(saveFlag,'Yes')
    
    %...if we don't have a file path
    if strcmpi(DPdata.file.path,'')
            
        %...get the current path and set it as DPdata.file.path...
        dummy=what;
        DPdata.file.path=dummy.path;
        
        %...and update the size of DPdata
        dummy=whos('DPdata');
        DPdata.file.size=dummy.bytes/1048576;%MBs=bytes/1024^2 = bytes/1048576 
        
    end
    
    %...bring together file name and path
    filePathName=[DPdata.file.path,'/',DPdata.file.name];
    %...if it alreadt exists...
    if exist(filePathName,'file')==2
        %...append it to the existing file
        save(filePathName, 'DPdata', '-append');
    else
        %...otherwise, just save it.
        save(filePathName, 'DPdata');
    end
    
end


%DONE!

function Axis=defineTimeFreqAxes(axisData,axisLabel,axisUnit,axisScale,axisValues)


...create unit...
axisUNIT = unit(axisUnit{1});

%...find minima and maxima and delta...
f = axisData(1);
if strcmpi(axisScale{1},'arbitrary')
    Axis = achse(axisUNIT,axisValues{1});
    Axis=setfirst(Axis,min(axisValues{1}));
    
else
    
    l = axisData(end);
    
    if strcmpi(axisScale{1},'linear')
        Nvals = length(axisData);
        d = (l-f)/(Nvals-1);
         %...create axis...
        Axis = achse(axisUNIT,f,d);
        
        lastlin=l;
    
    else
        
        if strcmpi(axisScale{1},'log2')
            %f = 2^floor(log2(axisData(1)));
            d=2;
             %l = round(log2(axisData(end)));
            %Nvals = 1+l-log2(f);
            
        elseif strcmpi(axisScale{1},'log10')
            %f = 10^floor(log10(axisData(1)));
            d=10;
            %l = round(log10(axisData(end)));
            %Nvals = 1+l-log10(f);
            
        elseif strcmpi(axisScale{1},'ln')
            %f = exp(floor(log(axisData(1))));
            d=exp(1);
            %l = round(log(axisData(end)));
            %Nvals = 1+l-log(f);
            
        end
        
        Nvals = ceil( (l-f)/d );
        firstlin = d^f;
        lastlin = d^l;
        %...create axis...
        Axis = achse(axisUNIT,firstlin,d,'log');
        
    end
    
    v=spacing(Axis,Nvals);
    
    LastIncluded = (v(end)>=lastlin);
    while ~LastIncluded
        Nvals=Nvals+1;
        v=spacing(Axis,Nvals);
        LastIncluded = (v(end)>=lastlin);
    end
    
    %...set values
    Axis = setvalues(Axis,v);
end
%...set its name...
Axis = setname(Axis,axisLabel{1});



function [dataAXES dataMINs dataMAXs dataDELTAs]=initializeDataAxes(axesLabels,axesUnits,axesScales,axesValues)

N=numel(axesLabels);

for jj=1:N;
    
    %...create unit...
    dataUNITS{jj} = unit(axesUnits{jj});
    
    dataAXES{jj} = achse(dataUNITS{jj});
    
    %...set scale...
    if strcmpi(axesScales{jj},'linear')
        dataAXES{jj} = setresolution(dataAXES{jj},'linear');
    elseif strcmpi(axesScales{jj},'arbitrary')
        dataAXES{jj} = achse(dataUNITS{jj},axesValues{jj});
        dataAXES{jj} = setfirst(dataUNITS{jj},axesValues{jj}(1));
    else
        dataAXES{jj} = setresolution(dataAXES{jj},'logarithmic');
    end
    %...set its name...
    dataAXES{jj} = setname(dataAXES{jj},axesLabels{jj});
    
    %We need the min and max of all datasets
    dataMINs(jj)=inf;
    dataMAXs(jj)=-inf;
    %...and the resolution as well for linear datasets
    if strcmpi(axesScales{jj},'linear')
        dataDELTAs(jj)=-inf;
    elseif strcmpi(axesScales{jj},'log2')
        dataDELTAs(jj)=2;
    elseif strcmpi(axesScales{jj},'log10')
        dataDELTAs(jj)=10;
    elseif strcmpi(axesScales{jj},'ln')
        dataDELTAs(jj)=exp(1);
    else
        dataDELTAs(jj)=[];
    end
    
end


function [Axes dataMINs dataMAXs dataDELTAs]=defineCurrDataAxes(dataAXES,dataMINs, dataMAXs, dataDELTAs,axesData,axesScales)

N=numel(dataAXES);

for jj=1:N;
    
    f=min(axesData{jj}(:));
    
    l=max(axesData{jj}(:));
    
    if f<dataMINs(jj)
        dataMINs(jj)=f;
    end
    if l>dataMAXs(jj)
        dataMAXs(jj)=l;
    end
    
    
    %Create a copy of dataAXIS for this dataset
    Axes{jj}=achse(dataAXES{jj});
    
    %...if the axis is not arbitrary and commong among all datasets...
    if ~strcmpi(axesScales{jj},'arbitrary')
        
        if strcmpi(axesScales{jj},'linear')
            
            diffs = diff(axesData{jj});
            d = max(diffs(:));
            if d>dataDELTAs(jj)
                dataDELTAs(jj)=d;
            end
            
            firstlin=f;
            lastlin=l;
            

        else
            
            d=dataDELTAs(jj);
            
%             if strcmpi(axesScales{jj},'log2')
%                 %f = 2^floor(log2(currMIN));
%                 d=2;
% %                 l = round(log2(currMAX));
% %                 Nvals = 1+l-log2(f);
%                 
%             elseif strcmpi(axesScales{jj},'log10')
%                 %f = 10^floor(log10(currMIN));
%                 d=10;
% %                 l = round(log10(currMAX));
% %                 Nvals = 1+l-log10(f);
%                 
%             elseif strcmpi(axesScales{jj},'ln')
%                 %f = exp(floor(log(currMIN)));
%                 d=exp(1);
% %                 l = round(log(currMAX));
% %                 Nvals = 1+l-log(f);
%
%             end

            firstlin = d^f;
            lastlin = d^l;

        end
        
        Nvals=ceil( (l-f)/d );
        Axes{jj} = setfirst(Axes{jj},firstlin);
        Axes{jj} = setdelta(Axes{jj},d);
        
        %...set values
        v=spacing(Axes{jj},Nvals);
        LastIncluded = (v(end)>=lastlin);

        while ~LastIncluded
            Nvals=Nvals+1;
            v=spacing(Axes{jj},Nvals);
            LastIncluded = (v(end)>=lastlin);
        end
        %...set values
        Axes{jj} = setvalues(Axes{jj},v);
        
    end
end


function dataAXES=updateDataAxes(dataAXES, dataMINs, dataMAXs, dataDELTAs,axesScales)

N=numel(dataAXES);

for jj=1:N;
    
  
    
    %Update also the general dataAXIS
    if ~strcmpi(axesScales{jj},'arbitrary')
        
        f = dataMINs(jj);
        d = dataDELTAs(jj);
        l = dataMAXs(jj);
        
        if strcmpi(axesScales{jj},'linear')
            %             f = dataMINs(jj);
            %             d = dataDELTAs(jj);
            %             l=dataMAXs(jj);
            %             Nvals=(l-f)/d + 1;
           firstlin=f;
           lastlin=l;
           
        else
            %             if strcmpi(axesScales{jj},'log2')
%             f = 2^floor(log2(dataMINs(jj)));
%             d=2;
%             l = round(log2(dataMAXs(jj)));
%             Nvals = 1+l-log2(f);
%             
%         elseif strcmpi(axesScales{jj},'log10')
%             f = 10^floor(log10(dataMINs(jj)));
%             d=10;
%             l = round(log10(dataMAXs(jj)));
%             Nvals = 1+l-log10(f);
%             
%         elseif strcmpi(axesScales{jj},'ln')
%             f = exp(floor(log(dataMINs(jj))));
%             d=exp(1);
%             l = round(log(dataMAXs(jj)));
%             Nvals = 1+l-log(f);
            firstlin = d^f;
            lastlin = d^l;
            
        end
        
        Nvals=ceil( (l-f)/d );
        
        dataAXES{jj} = setfirst(dataAXES{jj},firstlin);
        dataAXES{jj} = setdelta(dataAXES{jj},d);

        %...set values
        v=spacing(dataAXES{jj},Nvals);
        LastIncluded = (v(end)>=lastlin);
        while ~LastIncluded    
            Nvals=Nvals+1;
            v=spacing(dataAXES{jj},Nvals);
            LastIncluded = (v(end)>=lastlin);
        end
        %...set values
        dataAXES{jj} = setvalues(dataAXES{jj},v);
        
    end
end
