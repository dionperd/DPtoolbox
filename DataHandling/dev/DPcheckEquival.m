function [Total Result Info] = DPcheckEquival(DPdata,checks)

%This function checks the equivalence of the DPdata structures in the
%DPdata structure array in terms of the dimensions determined in the chechs
%argument

%Inputs:
%   -DPdata: DPdata structure arrat
%   -checks: cell string whose elements can take any of the following
%           values: 'domain','datasetLabels','fsample','time','freq','dims','tau','Ntrials',
%           'Nchannels','dataSIZE', 'axisLabels', 'axisUnits',
%           'axisSales','axisValues'

%Output:
%   -Total: 1 if all checks are true, 0 if at least one is false, -1 in
%   all other cases
%
%   Optional outputs:
%   -Result: structure containing fields indexed by the name of each check
%   and where Result.(checks{ii})= 1 for TRUE, 0 for FALSE
%   -Info: structure containing the elements of the DPdata structures that
%   were checked
%

Total=-1;
Result=struct();
Info=struct();

N_datasets = numel(DPdata);
if N_datasets<2
    warning('DPcheckEquival was terminated because DPdata containes less than 2 datasets')
    return
end

%Validate the inputs
%[RESULT x] = DPvalidateData(x,testfun,param,mode,execfun,default,varName,funcName)
funcName='DPcheckEquival';
varName='checks';
testCHECKS = { @(checks)iscellstr(checks),...
               @(checks) any( ( numel(checks)==[1:16] ) ),...
               @(checks) cellfun(@(check_i) any( strcmpi(checks,{...
                 'domain','fsample','Ntimes','time',...
                 'Nchannels','Ntrials','dataSIZE', 'dims',...
                 'datasetLabels','Nfreqs','freq','tau', ...
                 'axisLabels', 'axisUnits', 'axisSales','axisValues'}...
                  ) ),checks,'uniformoutput',false) };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[RESULT checks] = DPvalidateData(checks,testCHECKS,param,mode,execfun,default,varName,funcName);




%Check domain
if any(strcmpi('domain',checks))                       %domain=
  [Result.('domain'), Info.('domain')]=DPcheckDomain(arrayfun(@(dataset) dataset.domain, DPdata, 'uniformoutput',false));
end

%Check fsample
if any(strcmpi('fsample',checks))                       %fsample=
  [Result.('fsample'), Info.('fsample')]=DPcheckFsample(cell2mat(arrayfun(@(dataset) dataset.fsample, DPdata, 'uniformoutput',false)));
end

%Check Ntimes
time=[];
if any(strcmpi('Ntimes',checks))
    time = arrayfun(@(dataset) dataset.time, DPdata, 'uniformoutput',false);
    [Result.('Ntimes'), Info.('Ntimes')]=DPcheckNtimes(cellfun(@(time)length(time),time));
end

%Check time
if any(strcmpi('time',checks))
    
    if isempty(time)
        time = arrayfun(@(dataset) dataset.time, DPdata, 'uniformoutput',false);
        [ResultNtimes, Ntimes]=DPcheckNtimes(cellfun(@(time)length(time),time));
        clear Ntimes;
    else
        ResultNtimes=Result.('Ntimes');
    end
    
    if ResultNtimes
        [Result.('time'), Info.('time')]=DPcheckTime(cell2mat(time));
    else
        Result.('time') = 0;
        Info.('time') = time;
    end
    
    clear ResultNtimes;
end

%Check dataSIZE
if any(strcmpi('dataSIZE',checks))                                 %dataSIZE=
  [Result.('dataSIZE'), Info.('dataSIZE')]=DPcheckDataSIZE(arrayfun(@(dataset) size(dataset.signals), DPdata, 'uniformoutput',false));
end

%Check Nchannels
if any(strcmpi('Nchannels',checks))     
    if isfield(Result.('dataSIZE'))
        if Result.('dataSIZE')
            Result.('Nchannels')=1;
            Info.('Nchannels')=Info.('dataSIZE')(1);
        else
            [Result.('Nchannels'), Info.('Nchannels')]=DPcheckNchannels(Info.('dataSIZE')(:,1));
        end
    else
                                                                            %Nchannels=
        [Result.('Nchannels'), Info.('Nchannels')]=DPcheckNchannels(cell2mat(arrayfun(@(dataset) size(dataset.signals,1), DPdata, 'uniformoutput',false)));
    end
end


%Check Ntrials
if any(strcmpi('Ntrials',checks))      
    if isfield(Result.('dataSIZE'))
        if Result.('dataSIZE')
            Result.('Ntrials')=1;
            Info.('Ntrials')=Info.('dataSIZE')(2);
        else
            [Result.('Ntrials'), Info.('Ntrials')]=DPcheckNtrials(Info.('dataSIZE')(:,2));
        end
    else                                                            %Ntrials=
        [Result.('Ntrials'), Info.('Ntrials')]=DPcheckNtrials(cell2mat(arrayfun(@(dataset) size(dataset.signals,2), DPdata, 'uniformoutput',false)));
    end
end

%Check dims
if any(strcmpi('dims',checks))                              
    if isfield(Result.('dataSIZE'))
        ResultdataSIZE=Result.('dataSIZE');
    else
        [ResultdataSIZE, dataSIZE]=DPcheckDataSIZE(arrayfun(@(dataset) size(dataset.signals), DPdata, 'uniformoutput',false));
        clear dataSIZE;
    end
    
    if ResultdataSIZE                                               %dims=
        [Result.('dims'), Info.('dims')]=DPcheckDims(arrayfun(@(dataset) dataset.dims, DPdata, 'uniformoutput',false));
    else
        Result.('dims')=0;
        Info.('dims') = arrayfun(@(dataset) dataset.dims, DPdata, 'uniformoutput',false);
    end
    clear ResultdataSIZE;
end




%Check datasetLabels

        datasetLabels = reshape(arrayfun(@(dataset) dataset.datasetLabels, DPdata, 'uniformoutput',false),N_datasets,1);
        if all(strcmpi(datasetLabels{1,1},datasetLabels{2:3,1})) && all(strcmpi(datasetLabels{1,1},datasetLabels{2:3,1}))
            domainRes=1;
        end


    
Total = all(Total>0);


function [Result, domain]=DPcheckDomain(domain)

Result=all( strcmpi( domain{1},{domain{2:end}} ) );
if Result
    domain = domain{1};
end



function [Result, fsample]=DPcheckFsample(fsample)
Result= all(fsample(2:end)==fsample(1));
if Result
    fsample = fsample(1);
end


function [Result, Ntimes]=DPcheckNtimes(Ntimes)
Result= all(Ntimes(2:end)==Ntimes(1));
if Result
    Ntimes = Ntimes(1);
end

function [Result, time]=DPcheckTime(time)
N_datasets = size(time,2);
time1 = repmat(time(:,1),Ntimes,N_datasets-1);    
Result= all(time(:)==time1(:));
if Result
    time = time(:,1);
end

function [Result, dataSIZE]=DPcheckDataSIZE(dataSIZE)
dataSIZE = cell2mat(cellfun(@(dataSIZE_i)dataSIZE_i(1:2),dataSIZE,'uniformoutput',false));
dataSIZE = [dataSIZE(1:2:end-1).' dataSIZE(2:2:end).'];
Result= all(dataSIZE(:,1)==dataSIZE(1,1)) && all(dataSIZE(:,2)==dataSIZE(1,2));
if Result
    dataSIZE = dataSIZE(1,:);
end

function [Result, Nchannels]=DPcheckNchannels(Nchannels)
Result= all(Nchannels(2:end)==Nchannels(1));
if Result
    Nchannels = Nchannels(1);
end

function [Result, Ntrials]=DPcheckNtrials(Ntrials)
Result= all(Ntrials(2:end)==Ntrials(1));
if Result
    Ntrials = Ntrials(1);
end

function [Result, dims]=DPcheckDims(dims)
N_datasets = numel(dims);
dimsMat = cell2mat(dims(2:N_datasets-1));
dims1 = repmat(dims{1},1,N_datasets-1);
Result= all(dimsMat(:)==dims1(:));
if Result
    dims = dims{1};
end
