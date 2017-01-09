function [DPdata EEGLABinfo]=DPeeglab2DPdata(EEG,domain,chansORcomps,varargin)%fileName,filePath,saveFlag

%This function converts EEGLAB data structures into a time series or time-frequency DPdata structure
%For time-frequency data EEG has to include a field EEG.freqs with the
%frequency vector of the time - frequency decomposition
%and a field EEG.tfdata(channels x freqs x times x trials) for chansORcomps='channels' or EEG.icaactTF(channels x freqs x times x trials) for
%chansORcomps: 'components'
%
%Input:
%   -EEG: an array of EEGLAB EEG data structures 
%   -domain: 'time' or 'time-freq'
%   -chansORcomps: 'channels' or 'components'
%   Optional inputs:
%   -datasetLabels: cell array of strings of numel equal to numel(size(data)). 
%        Each string corresponds to the label of each index of DPdata.data{...}
%        Default is {'','',...,''}
%   -fileName: a string of the file name (with extension) where the DPdata should be saved
%   -filePath: a string of the file path
%   -saveFlag: "Yes" or "No" to save the structure or not

%
%Output: 
%   -DPdata: an array of DPdata structures
%   -EEGLABinfo: a structure array that contains all
%                EEG structures but with EEG.data=[], EEG.times=[],
%                (EEG.freqs=[]), EEG.icaact=[], (EEG.icaactTF=[])  



%Validate the inputs
%[RESULT x] = DPvalidateData(x,testfun,param,mode,execfun,default,varName,funcName)
funcName='DPeeglab2DPdata';

%Required inputs
varName='domain';
testDOMAIN = { @(domain)ischar(domain),...
               @(domain)isvector(domain),... 
               @(domain)any( strcmpi(domain,'time')||strcmpi(domain,'time-freq') ) };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default=nan;
[RESULT domain] = DPvalidateData(domain,testDOMAIN,param,mode,execfun,default,varName,funcName);


%Validate chansORcomps
varName='chansORcomps';
testCHANSorCOMPS = { @(chansORcomps)ischar(chansORcomps),...
                     @(chansORcomps)isvector(chansORcomps),... 
                     @(chansORcomps)strcmpi(chansORcomps,'channels')||strcmpi(chansORcomps,'components') };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default=nan;
[RESULT chansORcomps] = DPvalidateData(chansORcomps,testCHANSorCOMPS,param,mode,execfun,default,varName,funcName);

    
%Validate EEG
varName='EEG';
testEEG = { @(EEG) isstruct(EEG),...
            @(EEG) isvector(EEG)};
param={{},{}};
mode=['e','e'];
execfun={{},{}};
default=nan;
[RESULT EEG] = DPvalidateData(EEG,testEEG,param,mode,execfun,default,varName,funcName);

N_EEG = numel(EEG);
    
testEEG_i = { @(EEG_i) ~isempty(EEG_i),...
              @(EEG_i) isfield(EEG_i,'times'),...
              @(EEG_i) isfield(EEG_i,'srate'),...
              @(EEG_i) isfield(EEG_i,'chanlocs') };
param={{},{},{},{}};
mode=['e','e','e','e'];
execfun={{},{},{},{}};
default=nan;              

if strcmpi(domain,'time')
    if strcmpi(chansORcomps,'channels')
        
        testEEG_i(5) = { @(EEG_i) isfield(EEG_i,'data')};
 
    else
        testEEG_i(5) = { @(EEG_i) isfield(EEG_i,'icaact')};
   
    end
    
    param(5)={{}};
    mode(5)=['e'];
    execfun(5)={{}};
    
else
    if strcmpi(chansORcomps,'channels')
        
        testEEG_i(5:6) = { @(EEG_i) isfield(EEG_i,'freqs'),...
                           @(EEG_i) isfield(EEG_i,'tfdata')};
 
    else
        testEEG_i(5:6) = { @(EEG_i) isfield(EEG_i,'freqs'),...
                           @(EEG_i) isfield(EEG_i,'icaactTF') };
   
    end
    
    param(5:6)={{},{}};
    mode(5:6)=['e','e'];
    execfun(5:6)={{},{}};
end

for ii=1:N_EEG;
    
    varName=['EEG(',num2str(ii),')'];
    [RESULT EEG(ii)] = DPvalidateData(EEG(ii),testEEG_i,param,mode,execfun,default,varName,funcName);
end



if nargin>3
    fileName=varargin{1};
    
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


if nargin>4
    filePath=varargin{2};
    
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

if nargin>5
    saveFlag=varargin{3};
    
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

%Set dataset labels
datasetLabels(1:2) = {chansORcomps,'trials'};

%Set axes information
if strcmpi(domain,'time')
    axesLabels={'time', 'EEG'};
    axesUnits={'ms','uV'};
    axesScales={'linear', 'linear'};
else
    axesLabels={'time', 'frequency','real','imaginative'};
    axesUnits={'ms','Hz','',''};
    axesScales={'linear', 'linear', 'linear', 'linear'};
end


%Start the main job....
for ii=1:N_EEG;
    
    currEEG = EEG(1);
    
    %Delete the first EEG dataset to reduce memory requirements
    EEG(1)=[];
    
    %Unload EEG structure
    
    if strcmpi(domain,'time')
        if strcmpi(chansORcomps,'channels')
            EEGdata = currEEG.data;
        else
            EEGdata = currEEG.icaact;
            
        end
        
        N_signals = size(EEGdata,1);
        N_trials = size(EEGdata,3);
        N_times = size(EEGdata,2);
        
        signals = reshape( mat2cell( double(EEGdata), ones(1,N_signals), N_times, ones(1,N_trials)  ), [N_signals,N_trials] );
        
        signals = cellfun(@(signal_i)signal_i.',signals,'uniformoutput',false); 
        
    else
                
        if strcmpi(chansORcomps,'channels')
            EEGdata = currEEG.tfdata;
        else
            EEGdata = currEEG.icaactTF;
            
        end
        
        N_signals = size(EEGdata,1);
        N_trials = size(EEGdata,4);
        N_times = size(EEGdata,2);
        N_freqs = size(EEGdata,3);
        
        signals = reshape( mat2cell( double(EEGdata), ones(1,N_signals), N_times, N_freqs, ones(1,N_trials)  ), [N_signals,N_trials] );
        
        signals = cellfun(@(signal_i)reshape(signal_i,N_times,N_freqs),signals,'uniformoutput',false); 
     
    end
    
    
    if isfield(currEEG,'data')
        currEEG.data=[];
    end
    if isfield(currEEG,'icaact')
        currEEG.icaact=[];
    end
    if isfield(currEEG,'tfdata')
        currEEG.tfdata=[];
    end
    if isfield(currEEG,'icaactTF')
        currEEG.icaactTF=[];
    end
    
    
    %...and EEGLAB info structure
    EEGLABinfo(ii)=currEEG;
    

    %Create DPdata
    if strcmpi(domain,'time')
        DPdata(ii) = DPcreateData(domain,currEEG.srate,currEEG.times.',signals,datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales);%(tau/freq),datasetLabels,fileName,filePath,axesLabels,axesUnits,axesScales,axesValues
    else
        DPdata(ii) = DPcreateData(domain,currEEG.srate,currEEG.times.',signals,currEEG.freqs,datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales);%(tau/freq),datasetLabels,fileName,filePath,axesLabels,axesUnits,axesScales,axesValues
    end
    %DONE!!!
    
    clear EEGdata;
    
end

if length(fileName)>4
    fileName=[fileName(1:end-4),'_tm',fileName(end-3:end)];
else
    fileName='';
end
%....check if the file path exists...
if exist(filePath,'dir')~=7
    filePath='';
end

%Saving DPdata in the hard disk
%If we have a file name, lets save the DPdata at the hard disk
if strcmpi(saveFlag,'Yes')
    
    %...if we don't have a file path
    if strcmpi(filePath,'')
            
        %...get the current path and set it as DPdata.file.path...
        dummy=what;
        filePath=dummy.path;
        
    end
    
    %...bring together file name and path
    filePathName=[filePath,'/',fileName];
    %...if it alreadt exists...
    if exist(filePathName,'file')==2
        %...append it to the existing file
        save(filePathName, 'DPdata', '-append');
    else
        %...otherwise, just save it.
        save(filePathName, 'DPdata');
    end
    
end

