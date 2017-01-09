function MergedData=DPmerge(DPdata,subset,output)

%This function merges different DPdata structures from a DPdata structure
%array either to a single DPdata structure or to a single cell of tstool signals or 
% to a single cell of data arrays, accross channels or trials. 
% The respective domains, sizes and axes have to be equivalent.
%DPdata can result in the output only if domains are "time"

%Input: 
%   -DPdata: DPdata structure array
%   -subset: a string that takes values 'channels' or 'trials' or
%            determining the subset of the datasets accross which they will be merged
%   -output: "DPdata', tstool' or 'array'
  
%Output: 
%   -MergedData: output DPdata structure or cell with tstool signals or
%                with data arrays



%Input data validation
funcName='DPmerge';

varName='subset';
testSUBSET = {@(subset)ischar(subset),...
              @(subset)isvector(subset),...
              @(subset) strcmpi(subset,'channels')||strcmpi(subset,'trials') };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[RESULT subset] = DPvalidateData(subset,testSUBSET,param,mode,execfun,default,varName,funcName);

varName='output';
testOUTPUT = {@(output)ischar(output),...
              @(output)isvector(output),...
              @(output) strcmpi(output,'DPdata')||strcmpi(output,'tstool')||strcmpi(output,'array') };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[RESULT output] = DPvalidateData(output,testOUTPUT,param,mode,execfun,default,varName,funcName);


%Check if all datasets are of the same domain
domain = arrayfun(@(dataset) dataset.domain, DPdata, 'uniformoutput',false);
if ~all(strcmpi(domain,domain{1}))
    error('Datasets cannot be merged because they are not of the same domain');
else
    domain=domain{1};
end
if strcmpi(domain,'time-freq')
    Naxes=4;
else
    Naxes=2;
end

%Check if time has the same length in all datasets
time = arrayfun(@(dataset) dataset.time, DPdata, 'uniformoutput',false);
Ntimes=cellfun(@(time)length(time),time);
if all(Ntimes==Ntimes(1))
    Ntimes=Ntimes(1);
else
    error('Datasets cannot be merged because they do not have the length in time');
end

%Check if the domain is 'time' to get DPdata in the output
if strcpmi(output,'DPdata')
    if ~strcpmi(domain,'time')
        error('Domain has to be ''time'' to get a DPdata structure in the output')
    end
end
        
%Check if all datasets have the correct sizes
N_datasets=numel(DPdata);
sizes = cell2mat(reshape(arrayfun(@(dataset) size(dataset.signals), DPdata, 'uniformoutput',false),N_datasets,1));
if strcmpi(subset,'channels')
    N_trials=sizes(1,2); %The number of trials has to be the same in order to merge channels
    if all(sizes(:,2)==N_trials)
        N_channels = sum(sizes(:,1)) %Number of channels if we merge them
    else
        error('Datasets cannot be merged because they do not have the same number of trials');
    end  
end
if strcmpi(subset,'trials')
    N_channels = sizes(1,1); %The number of channels has to be the same in order to merge trials
    if all(sizes(:,1)==N_channels)
        N_trials = sum(sizes(:,2)) %Number of trials if we merge them
    else
        error('Datasets cannot be merged because they do not have the same number of channels');
    end
end

if strcpmi(output,'DPdata')

%Check if time and axes are equivalent if we want to get a DPdata in the
%output
time = cellmat(time);
time1 = repmat(time(:,1),Ntimes,N_datasets);
if ~all(time(:)==time1(:))
    time=time(:,1);
    clear time1;
else
    error('time has to be the same accross datasets to get a DPdata structure in the output')
end

[axesLabels1, axesUnits1, axesScales1, axesValues1] = DPgetAxisInfo(DPdata(1).axes);

[axesLabels, axesUnits, axesScales, axesValues] = arrayfun(@(axCell)DPgetAxisInfo(axCell),DPdata(2:end).axes, 'uniformoutput',false);


        
    

