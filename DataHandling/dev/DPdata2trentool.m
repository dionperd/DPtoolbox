function trentoolDATA=DPdata2trentool(DPdata,channelLabels)

%This function creates trentool data structures 
%from the datasets with indexes in the range of Ranges of a time series DPdata
%structure
%
%Input:
%   -DPdata: the DPdata structure array
%   -channelLabels: cell string with the labels of the channels/components
%   etc
%
%Output: 
%   -trentoolDATA: a trentool structure array:

%

%Validate ind:
funcName='DPdata2trentool';

varName='channelLabels';
testLABELS = { @(channelLabels)iscellstr(channelLabels),...
               @(channelLabels)isvector(channelLabels)     };
param={{},{}};
mode=['e','e'];
execfun={{},{}};
default=nan;
[RESULT channelLabels] = DPvalidateData(channelLabels,testLABELS,param,mode,execfun,default,varName,funcName);
N_channelLabels = numel(channelLabels);

testLABELS_i= { @(channelLabel_i)ischar(channelLabel_i),...
                @(channelLabel_i)isvector(channelLabel_i),...
                @(channelLabel_i)~isempty(channelLabel_i)};
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default=nan;
for ii=1:N_channelLabels;
    varName=['channelLabels{',num2str(ii),'}'];
    [RESULT channelLabels{ii}] = DPvalidateData(channelLabels{ii},testLABELS_i,param,mode,execfun,default,varName,funcName);
end


for ii=1:numel(DPdata);
    
    dataSIZE=size(DPdata(1).signals);
    N_Channels = dataSIZE(1);
    if N_channelLabels~=N_Channels;
        error('N_Channels of DPdata(%d) is not equal to N_channelLabels',ii)
    end
    N_Trials = dataSIZE(2);
    

    %All datasets have to be of dimension 1
    varName=['DPdata(',num2str(ii),'.dims'];
    testDIMS= { @(dims) all( arrayfun(@(dim)any(dim==[0 1]),dims(:)) ) };
    param={{}};
    mode=['e'];
    execfun={{}};
    default=nan;
    [RESULT DPdata(1).dims] = DPvalidateData(DPdata(1).dims,testDIMS,param,mode,execfun,default,varName,funcName);
    
    Ntimes=length(DPdata(1).time);
    
    %Move the selected DPdata.data into a trentool structure data
    for jj=1:N_Trials;
        for kk=1:N_Channels
            trentoolDATA(ii).trial{jj}(kk,:) = data( DPdata(1).signals{kk,jj} );
        end
    end
    
    %Add time
    trentoolDATA(ii).time = mat2cell([repmat(DPdata(1).time,1,N_Trials)],Ntimes,ones(1,N_Trials));
    
    %Add fsample
    trentoolDATA(ii).fsample = DPdata(1).fsample;
    
    %Add channel labels
    trentoolDATA(ii).label = channelLabels;
    
    
    
    DPdata(1)=[];
    
end



