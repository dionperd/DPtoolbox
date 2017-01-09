function linInd = DPsubInd2linIndRange(siz,Ranges)

%This function converts a range of subindexes to linear indexes
%
%Inputs: 
%   -siz: vector that corresponds to the size of the indexed array
%   -Ranges: cell of vectors that determins the ranges of values per each axis
%
%Outputs:
%   -linInd: linear indexes sorted in ascending order

%Validate ind:
funcName='DPsubInd2linIndRange';

varName='siz';
testSIZ = { @(siz)isnumeric(siz),...
            @(siz)isvector(siz),...
            @(siz)(~isempty(siz)),...
            @(siz)isreal(siz),...
            @(siz)all(siz>0),...
            @(siz)all(round(siz)==siz)};
param={{},{},{},{},{},{}};
mode=['e','e','e','e','e','e'];
execfun={{},{},{},{},{},{}};
default=nan;
[RESULT siz] = DPvalidateData(siz,testSIZ,param,mode,execfun,default,varName,funcName);
Nsiz=length(siz);


varName='Ranges';
testRANGES = { @(Ranges)iscell(Ranges),...
              @(Ranges,Nsiz)(numel(Ranges)==Nsiz) };
param={{},{Nsiz}};
mode=['e','e'];
execfun={{},{}};
default=nan;
[RESULT Ranges] = DPvalidateData(Ranges,testRANGES,param,mode,execfun,default,varName,funcName);

for ii=1:Nsiz;
    
    varName=['Ranges{',num2str(ii),'}'];
    testRANGE_i = { @(Range_i)isnumeric(Range_i),...
                    @(Range_i)isvector(Range_i),...
                    @(Range_i)(~isempty(Range_i)),...
                    @(Range_i)( length(Range_i)==length(unique(Range_i)) ),...
                    @(Range_i,siz)( all( arrayfun( @(ind)any(ind==1:siz), Range_i ) ) ) };
    param={{},{},{},{},{siz(ii)}};
    mode=['e','e','e','e','e'];
    execfun={{},{},{},{},{}};
    default=nan;
    [RESULT Ranges{ii}] = DPvalidateData(Ranges{ii},testRANGE_i,param,mode,execfun,default,varName,funcName);
    
end


switch Nsiz
    
    case 1 
        
        linInd = sort(Ranges{1});
        
    %linearInd = sub2ind(arraySize, dim1Sub, dim2Sub, dim3Sub, ...)
    case 2 
        
        %[X,Y] = meshgrid(x,y)
        [subInd{1} subInd{2}] =meshgrid(Ranges{1},Ranges{2});
        for ii=1:2;
            subInd{ii}=subInd{ii}(:);
        end
        
        linInd = sort(sub2ind(siz,subInd{1}, subInd{2}));
        
    case 3
        
        %[X,Y,Z] = meshgrid(x,y,z)
        [subInd{1} subInd{2} subInd{3}] =meshgrid(Ranges{1},Ranges{2},Ranges{3});
        for ii=1:3;
            subInd{ii}=subInd{ii}(:);
        end
        
        linInd = sort(sub2ind(siz,subInd{1}, subInd{2}, subInd{3}));
        
    otherwise
        
        %[X1,X2,X3,...] = ndgrid(x2,x1,x3,...)
        COMMAND=['[subInd{2}, subInd{1}'];
        for ii=3:Nsiz;
            COMMAND=[COMMAND,', subInd{',num2str(ii),'}'];
        end
        COMMAND=[COMMAND,']=ndgrid(Ranges{2}, Ranges{1}'];
        for ii=3:Nsiz;
            COMMAND=[COMMAND,', Ranges{',num2str(ii),'}'];
        end
        COMMAND=[COMMAND,');'];
        eval(COMMAND);
        
        for ii=1:Nsiz;
            COMMAND = ['subInd{',num2str(ii),'}=subInd{',num2str(ii),'}(:);'];
            eval(COMMAND);
        end
        COMMAND=['linInd=sort(sub2ind(siz'];
        for ii=1:Nsiz;
            COMMAND=[COMMAND,', subInd{',num2str(ii),'}'];
        end
        COMMAND=[COMMAND,'));'];
        eval(COMMAND);
end