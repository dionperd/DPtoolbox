function DPdata=DPscaleTransform(DPdata,axesScales,varargin)%filename, filepath,saveFlag

%This function transforms the DPdata axes' scales and data accordingly.
%
%Input: 
%   -DPdata: DPdata structure array of equivalent DPdata structures 
%   -axesScales: cell string determining the desired new scales for each axis
%                 it takes values ('linear','log2', 'log10' or 'ln'))
%   Optional inputs:
%   -filename
%   -filepath
%   -saveFlag
%
%Output: 
%   -DPdata: transformed DPdata structure 

%Input data validation
funcName='DPloglinTransform';

N_datasets = numel(DPdata);

%Initialization
%Defaults are the existing ones...
if strcmpi(DPdata(1).domain,'time-freq')
    Naxes=4;
else
    Naxes=2;
end
    
varName='axesScales';
testAXESSCALES = { @(axesScales)iscellstr(axesScales),...
                   @(axesScales,Naxes)( size(axesScales,2)==Naxes ),...
                   @(axesScales,N_datasets)any( ( size(axesScales,1)==[N_datasets,1] ) )};

param={{},{Naxes},{N_datasets}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[RESULT axesScales] = DPvalidateData(axesScales,testAXESSCALES,param,mode,execfun,default,varName,funcName);

NaxesScales = size(axesScales,1);
if RESULT==1
    
    for ii=1:NaxesScales;
        for jj=1:Naxes;
            
            varName=['axesScales{',num2str(ii),',',num2str(jj),'}'];
            testAXESSCALE_i = { @(axesScale_i)strcmpi(axesScale_i,'linear')||strcmpi(axesScale_i,'log2')||strcmpi(axesScale_i,'log10')||strcmpi(axesScale_i,'ln')||strcmpi(axesScale_i,'arbitrary')};
            param={{}};
            mode=['e'];
            execfun={{}};
            default='';
            [RESULT axesScales{ii,jj}] = DPvalidateData(axesScales{ii,jj},testAXESSCALE_i,param,mode,execfun,default,varName,funcName);
            
        end
    end
    
end
if (NaxesScales==1)
    axesScales=repmat(axesScales,N_datasets,1);
end

if nargin>2
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


if nargin>3
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

if nargin>4
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

                  %linear/        
                  %arbitrary   log2                    log10                  ln                      
transformFuncs = { @(x)x,      @(x)log2(x),            @(x)log10(x),          @(x)log(x);... 
                   @(x)2.^x,   @(x)x,                  @(x)log2(x)./log2(10), @(x)log2(x)./log2(exp(1));...
                   @(x)10.^x,  @(x)log10(x)./log10(2), @(x)x,                 @(x)log10(x)./log10(exp(1));...
                   @(x)exp(x), @(x)log(x)./log(2),     @(x)log(x)./log(10),   @(x)x                          };

%Find the index of the scale given that arbitrary and linear are considered
%both to be in a linear space
findScaleInd = {@(scale)any(strcmpi(scale,{'linear','arbitrary'})), @(scale)strcmpi(scale,'log2'), @(scale)strcmpi(scale,'log10'), @(scale)strcmpi(scale,'ln') }; 



for ii=1:N_datasets;
    
    [axesLabels, axesUnits, defaultAxesScales, axesValues] = DPgetAxisInfo(DPdata(ii).axes);

    %Time axis:
    TIMEscaleFROMind = find(cellfun(@(fun)fun(defaultAxesScales{1}),findScaleInd));
    TIMEscaleTOind = find(cellfun(@(fun)fun(axesScales{ii,1}),findScaleInd));
    
    if strcmpi(DPdata(ii).domain,'time-freq')
        %Frequency axis:
        FREQscaleFROMind = find(cellfun(@(fun)fun(defaultAxesScales{2}),findScaleInd));
        FREQscaleTOind = find(cellfun(@(fun)fun(axesScales{ii,2}),findScaleInd));
        for jj=3:4;
            %If we need to transform the data...
            DATAscaleFROMind(jj-2) = find(cellfun(@(fun)fun(defaultAxesScales{jj}),findScaleInd));
            DATAscaleTOind(jj-2) = find(cellfun(@(fun)fun(axesScales{ii,jj}),findScaleInd));
        end
    else
        %Data axis:
        DATAscaleFROMind = find(cellfun(@(fun)fun(defaultAxesScales{2}),findScaleInd));
        DATAscaleTOind = find(cellfun(@(fun)fun(axesScales{ii,2}),findScaleInd));
    end

    %If we need to transform time...
    if TIMEscaleFROMind~=TIMEscaleTOind
        %do so...
        DPdata(ii).time = transformFuncs{TIMEscaleFROMind,TIMEscaleTOind}(DPdata(ii).time);
    end
    
    %Get the data of the signals
    datasets=cellfun(@(signal_i) data(signal_i),DPdata(ii).signals,'UniformOutput',false);
    
    if strcmpi(DPdata(ii).domain,'time-freq')
        
        %If we need to transform frequency...
        if FREQscaleFROMind~=FREQscaleTOind
            %do so...
            DPdata(ii).freq = transformFuncs{FREQscaleFROMind,FREQscaleTOind}(DPdata(ii).freq);
            
        end
        
        %Data axes:
        for jj=3:4;
            if DATAscaleFROMind(jj-2)~=DATAscaleTOind(jj-2)
                %do so...
                datasets = cellfun(@(dataset) transformFuncs{DATAscaleFROMind(jj-2),DATAscaleTOind(jj-2)}(dataset(:,:,jj-2)), datasets,'UniformOutput',false);
            end
            
        end
    else        
        %If we need to transform the data...
        if DATAscaleFROMind~=DATAscaleTOind
            %do so...
            datasets = cellfun(@(dataset)transformFuncs{DATAscaleFROMind,DATAscaleTOind}(dataset),datasets,'UniformOutput',false);
        end
        
    end
    
    
    if strcmpi(DPdata(ii).domain,'time')
        DPdata(ii) = DPcreateData(DPdata(ii).domain,DPdata(ii).fsample,DPdata(ii).time,datasets,DPdata(ii).datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales(ii,:),axesValues);
    elseif strcmpi(DPdata(ii).domain,'time-freq')
        DPdata(ii) = DPcreateData(DPdata(ii).domain,DPdata(ii).fsample,DPdata(ii).time,datasets,DPdata(ii).freq,DPdata(ii).datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales(ii,:),axesValues);
    else
        DPdata(ii) = DPcreateData(DPdata(ii).domain,DPdata(ii).fsample,DPdata(ii).time,datasets,DPdata(ii).tau,DPdata(ii).datasetLabels,fileName,filePath,'No',axesLabels,axesUnits,axesScales(ii,:),axesValues);
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
