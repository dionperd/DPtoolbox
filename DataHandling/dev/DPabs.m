function DPdata=DPabs(DPdata,varargin)%filename,filepath,saveFlag

%This function replaces DPdata with its absolute values 
%
%Input: 
%   -DPdata: DPdata structure array

%   Optional inputs:
%   -filename
%   -filepath
%   -saveFlag
%Output: 
%   -DPdata: transformed DPdata structure array


%Input data validation
funcName='DPabs';

if nargin>1
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


if nargin>2
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

if nargin>3
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

for ii=1:numel(DPdata);
    
    %Initialization
    [axesLabels, axesUnits, axesScales, axesValues] = DPgetAxisInfo(DPdata(ii).axes);  
    
    
    %Get the abs of the data of the signals
    datasets=cellfun(@(signal_i) abs(data(signal_i)),DPdata(ii).signals,'UniformOutput',false);
    
    
    if strcmpi(DPdata(ii).domain,'time')
        DPdata(ii) = DPcreateData(DPdata(ii).domain,DPdata(ii).fsample,DPdata(ii).time,datasets,DPdata(ii).datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales,axesValues);
    elseif strcmpi(DPdata(ii).domain,'time-freq')
        DPdata(ii) = DPcreateData(DPdata(ii).domain,DPdata(ii).fsample,DPdata(ii).time,datasets,DPdata(ii).freq,DPdata(ii).datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales,axesValues);
    else
        DPdata(ii) = DPcreateData(DPdata(ii).domain,DPdata(ii).fsample,DPdata(ii).time,datasets,DPdata(ii).tau,DPdata(ii).datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales,axesValues);
    end

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