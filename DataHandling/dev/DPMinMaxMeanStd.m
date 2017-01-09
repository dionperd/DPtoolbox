function Results=DPMinMaxMeanStd(DPdata,measures,subset)

%This function calculates any of the min, max, mean and std values accross the subset
%of DPdata structure array defined by the argument subset
%
%Input: 
%   -DPdata: DPdata structure array
%   -subset: a string that takes values 'datasets', 'signals', 'channels', 'trials' or
%                    determining the subset of the
%                   datasets accross which in, max, mean and variance should be calculated
%   -DatasetConcad: "Yes" or "No" for concadenating datasets of DPdata
%                    structure array or not, default='No'
%   -measures: any or all from 'mean','std','min','max'
  
%Output: 
%   -Results: output structure



%Input data validation
funcName='DPMinMaxMeanStd';

varName='measures';
testMEASURES = {@(measures)iscellstr(measures),...
                @(measures)isvector(measures),...   
                @(measures)all(cellfun(@(measure_i) ischar(measure_i)&&isvector(measure_i), measures),'uniformoutput',false)   };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[RESULT measures] = DPvalidateData(measures,testMEASURES,param,mode,execfun,default,varName,funcName);

funcs=struct();
Nfuncs=0;
funcNAMES = {};
if any(strcmpi('mean',measures))
    Nfuncs=Nfuncs+1;
    funcNAMES(Nfuncs)='mean';
    funcs.('mean') = @(x)mean(x(:));
end
if any(strcmpi('std',measures))
    Nfuncs=Nfuncs+1;
    funcNAMES(Nfuncs)='std';
    funcs.('std') = @(x)std(x(:));
end
if any(strcmpi('min',measures))
    Nfuncs=Nfuncs+1;
    funcNAMES(Nfuncs)='min';
    funcs.('min') = @(x)min(x(:));
end
if any(strcmpi('max',measures))
    Nfuncs=Nfuncs+1;
    funcNAMES(Nfuncs)='max';
    funcs.('max') = @(x)max(x(:));
end
 
if isempty(funcs)
    error('None of the mean, std, min or max measures was selected!')
end

varName='subset';
testSUBSET = {@(subset)ischar(subset),...
              @(subset)isvector(subset),...
              @(subset) strcmpi(subset,'datasets')||strcmpi(subset,'channels')||strcmpi(subset,'trials')||strcmpi(subset,'signals') };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[RESULT subset] = DPvalidateData(subset,testSUBSET,param,mode,execfun,default,varName,funcName);

varName='DatasetConcad';
testSUBSET = {@(DatasetConcad)ischar(DatasetConcad),...
              @(DatasetConcad)isvector(DatasetConcad),...
              @(DatasetConcad)strcmpi(DatasetConcad,'Yes')||strcmpi(DatasetConcad,'No') };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='No';
[RESULT subset] = DPvalidateData(subset,testSUBSET,param,mode,execfun,default,varName,funcName);


N_datasets=numel(DPdata);
dataSIZEs = reshape(cell2mat(arrayfun(@(dataset) size(dataset.signals), DPdata, 'uniformoutput',false)),N_datasets,2);

if strcmpi(subset,'DatasetConcad')
    
    if strcmpi(subset,'signals')
        
    
    Results = cell(1,Nfuncs);
    
    Results
    
else
    Results = cell(N_datasets,Nfuncs);
    
    
    
    
end




for ii=1:N_datasets;
    
   

end

