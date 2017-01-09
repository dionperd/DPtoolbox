function DPdata=DPnormal(DPdata,normalmode,cfg,varargin)%filename,filepath,saveFlag

%This function replaces DPdata with its absolute values 
%
%Input: 
%   -DPdata: DPdata structure array
%   -normalmode: 'center', 'normal', 'linear'
%   -cfg: a structure with fields
%       for 'linear':
%           -minVal, scalar, default =0 
%           -maxVal, scalar, default =1 
%           -subset: a string that takes values 'all','datasets', 'channels', 'trials' or
%                   'signals' determining the subset of the
%                   datasets within which min and max should be
%                   calculated
%       for 'normal': 
%           -mean
%           -std
%       or  -subset: a string that takes values 'all','datasets', 'channels', 'trials' or
%                   'signals' determining the subset of the
%                   datasets within which mean and variance should be calculated
%       for 'center': 
%           -mean
%       or  -subset: a string that takes values 'all','datasets', 'channels', 'trials' or
%                   'signals' determining the subset of the
%                   datasets within which mean should be
%                   calculated

%   Optional inputs:
%   -filename
%   -filepath
%   -saveFlag
%   
%Output: 
%   -DPdata: transformed DPdata structure array



%Input data validation
funcName='DPnormal';

varName='normalmode';
testNORMALMODE = {@(normalmode)ischar(normalmode),...
                  @(normalmode)isvector(normalmode),...
                  @(normalmode)strcmpi(normalmode,'center')||strcmpi(normalmode,'normal')||strcmpi(normalmode,'linear')};
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[RESULT normalmode] = DPvalidateData(normalmode,testNORMALMODE,param,mode,execfun,default,varName,funcName);

varName='cfg';

if strcmpi(normalmode,'linear')
    
    testCFG = {@(cfg)isstruct(cfg),...
               @(cfg)isfield(cfg,'subset'),...
               @(cfg)ischar(cfg.subset),...
               @(cfg)isvector(cfg.subset),...
               @(cfg)strcmpi(cfg.subset,'all')||strcmpi(cfg.subset,'datasets')||strcmpi(cfg.subset,'channels')||strcmpi(cfg.subset,'trials')||strcmpi(cfg.subset,'signals') };
    param={{},{},{},{},{}};
    mode=['e','e','e','e','e'];
    execfun={{},{},{},{},{}};
    default='';
    [RESULT cfg] = DPvalidateData(cfg,testCFG,param,mode,execfun,default,varName,funcName);
    
    if isfield(cfg.minVal) && isfield(cfg.maxVal) 
        varName='cfg.minVal & cfg.maxVal';
        testMINMAXVAL = {@(cfg)isnumeric(cfg.minVal),...
                         @(cfg)isnumeric(cfg.maxVal),...
                         @(cfg)isscalar(cfg.minVal),...
                         @(cfg)isscalar(cfg.maxVal),...
                         @(cfg)isreal(cfg.minVal),...
                         @(cfg)isreal(cfg.maxVal),...
                         @(cfg)cfg.maxVal>cfg.minVal};
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='';
        [RESULT cfg] = DPvalidateData(cfg,testMINMAXVAL,param,mode,execfun,default,varName,funcName);
    else
        cfg.minVal=0;
        cfg.maxVal=1;
    end

    
elseif strcmpi(normalmode,'normal')
    
    testCFG = {@(cfg)isstruct(cfg),...
               @(cfg)isfield(cfg,'subset')||(isfield(cfg,'mean')&&isfield(cfg,'std')) };
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    default='';
    [RESULT cfg] = DPvalidateData(cfg,testCFG,param,mode,execfun,default,varName,funcName);
    
    if isfield(cfg.subset) 
        varName='cfg.subset';
        testSUBSET = {@(cfg)ischar(cfg.subset),...
                      @(cfg)isvector(cfg.subset),...
                      @(cfg)strcmpi(cfg.subset,'all')||strcmpi(cfg.subset,'datasets')||strcmpi(cfg.subset,'channels')||strcmpi(cfg.subset,'trials')||strcmpi(cfg.subset,'signals') };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='';
        [RESULT cfg] = DPvalidateData(cfg,testSUBSET,param,mode,execfun,default,varName,funcName);
        
    elseif (isfield(cfg,'mean')&&isfield(cfg,'std'))
        varName='cfg.mean & cfg.std';
        testMEANSTD = {@(cfg)isnumeric(cfg.mean),...
                       @(cfg)isnumeric(cfg.std),...
                       @(cfg)isscalar(cfg.mean),...
                       @(cfg)isscalar(cfg.std),...
                       @(cfg)isreal(cfg.mean),...
                       @(cfg)isreal(cfg.std),...
                       @(cfg)cfg.std>0              };
        param={{},{},{},{},{},{},{}}; 
        mode=['e','e','e','e','e','e','e'];
        execfun={{},{},{},{},{},{},{}};
        default='';
        [RESULT cfg] = DPvalidateData(cfg,testMEANSTD,param,mode,execfun,default,varName,funcName);
                   
    end
    
   
    
elseif strcmpi(normalmode,'central')
    
    testCFG = {@(cfg)isstruct(cfg),...
               @(cfg)isfield(cfg,'subset')||isfield(cfg,'mean')};
    param={{},{}};
    mode=['e','e'];
    execfun={{},{}};
    default='';
    [RESULT cfg] = DPvalidateData(cfg,testCFG,param,mode,execfun,default,varName,funcName);
    
    if isfield(cfg.subset) 
        varName='cfg.subset';
        testSUBSET = {@(cfg)ischar(cfg.subset),...
                      @(cfg)isvector(cfg.subset),...
                      @(cfg)strcmpi(cfg.subset,'all')||strcmpi(cfg.subset,'datasets')||strcmpi(cfg.subset,'channels')||strcmpi(cfg.subset,'trials')||strcmpi(cfg.subset,'signals') };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='';
        [RESULT cfg] = DPvalidateData(cfg,testSUBSET,param,mode,execfun,default,varName,funcName);
        
    elseif isfield(cfg,'mean')
        varName='cfg.mean';
        testMEAN = {@(cfg)isnumeric(cfg.mean),...
                    @(cfg)isscalar(cfg.mean),...
                    @(cfg)isreal(cfg.mean)    };
        param={{},{},{}}; 
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='';
        [RESULT cfg] = DPvalidateData(cfg,testMEAN,param,mode,execfun,default,varName,funcName);
    end
    
end

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


N_datasets=numel(DPdata);

if strcmpi(normalmode,'linear')
        
    %Define the utility functions
    utilFunc = {@(x)min(x(:)), @(x)max(x(:))};
    
    %Define the normalization function
    normalFunc = @(x,minVal,maxVal,x_min,x_max)( (maxVal-minVal)./(x_max-x_min) ) .* (x-x_min) + minVal;
    
elseif strcmpi(normalmode,'normal')
    
     %Define the utility functions
    utilFunc = {@(x)mean(x(:)), @(x)std(x(:))};

    %Define the normalization function
    normalFunc = @(x,mean,std)(x-mean)./std;
    
elseif strcmpi(normalmode,'central')

    %Define the utility function
    utilFunc = @(x)mean(x(:));
    
    %Define the normalization function
    normalFunc = @(x,mean)x-mean;
    
end




for ii=1:N_datasets;
    
    %Initialization
    for jj=1:numel(DPdata(ii).axes);
        axesLabels{jj} = name(DPdata(ii).axes{jj});
        axesUnits{jj} = char(unit(DPdata(ii).axes{jj}));
        axesScales{jj} = resolution(DPdata(ii).axes{jj});
        axesValues{jj} = values(DPdata(ii).axes{jj});
    end

    
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