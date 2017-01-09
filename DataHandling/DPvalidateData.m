function [RESULT x] = DPvalidateData(x,testfun,param,mode,execfun,default,varName,funcName)

%This function:
%  1. checks whether argument x satisfies requirements determing by
%  functions testfun(x,param) with parameters param
%  2. discriminates between warnings and errors
%  3. has the possibility for correcting x in the case of warnings by
%  applying the correcting function execfun(x,param)
%  4. prints the respective colored messages
%  5. sustitutes x with default value or...
%  6. ...stops the execution with an error message in case of validation
%     failure, if default=nan
%
%  Inputs: 
%  -x, the argument to be validated
%  -testfun: a cell of function handles of functions to be tested 
%  -param: a cell of parameters (param{i} is empty if testfun{i} does need
%          any parameters), each param{i} may be a cell containing more
%          than one parameters
%  -mode: a vector of characters of value 'e' for errors,'w' for warnings
%         and 'c' for warnings with correction
%  -execFun: a cell of function handles of functions to be executed for the
%  correction of x in the case of warnings
%  -varName: a string determining the variable x name
%  -funcName: a string determining the function's name that validates x

%Create the input parser
p=inputParser;
p.FunctionName = 'DPvalidateData';
p.CaseSensitive=false; %NOT case sensitive
p.KeepUnmatched = false; %do not accept inputs undeclared here
p.StructExpand = false; %accept structures as single inputs

%Define the inputs and their classes and attributes:
%Required inputs
p.addRequired('x',@(x)1);
p.addRequired('testfun',@(testfun)iscell(testfun)&&all(cellfun(@(x)isa(x,'function_handle'),testfun)));
Ntests=numel(testfun);
p.addRequired('param',@(param)iscell(param)&&(numel(param)==Ntests));
p.addRequired('mode',@(mode)ischar(mode)&&(numel(mode)==Ntests));
p.addRequired('execfun',@(execfun)iscell(execfun) && all(cellfun(@(x)isa(x,'function_handle')||isempty(x),execfun)) && (numel(execfun)==Ntests));
p.addRequired('default',@(default)1);
p.addRequired('varName',@(varName)ischar(varName));
p.addRequired('funcName',@(funcName)ischar(funcName));

%Check inputs 
p.parse(x,testfun,param,mode,execfun,default,varName,funcName);


%For each test as long we don't get an error...
i=1;
CONTINUE=true;
while CONTINUE
    
    %...apply test... 
    if isempty(param{i})
        %...without...
        RESULT(i) = single(testfun{i}(x));
    else
        %...or with parameters...
        RESULT(i) = single(testfun{i}(x,param{i}{:}));

    end
    
    %...if the test fails...
    if (RESULT(i)==0)
       %...get the name of the function that failed...
        dummyf = functions(testfun{i});
        
        %...if it is an error...
        if strcmpi(mode(i),'e')
            
            %...signal it as such...
            RESULT(i)=-1;
            if isnumeric(default)&&all(isnan(default))
                %...print message...
                %cprintf('ERROR',['\nERROR in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\n'])
                fprintf(['\nERROR in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\n']);
                error('Function %s failed',funcName);
            else
                x=default;
                %cprintf('TEXT',['\nERROR in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\n'])
                fprintf(['\nERROR in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\n'])
                %cprintf('TEXT',['Function ',funcName,' continues setting default value for ',varName,'.\n'])
                fprintf(['Function ',funcName,' continues setting default value for ',varName,'.\n'])

            end
            
            
            %...if it is a warning...
        else
            %...signal it as such...
            RESULT(i)=0;
            
            %...and if x should be corrected...
            if ~isempty(execfun{i})
                %...get the name of the corection function...
                dummyf2 = functions(execfun{i});
                %...print message...
                %cprintf('TEXT',['\nWARNING in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\nExecuting correcting function ', dummyf2.function,'\n'])
                fprintf(['\nWARNING in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\nExecuting correcting function ', dummyf2.function,'\n'])                
                %...apply corection... 
                if isempty(param{i})
                    %...without...
                    x=execfun{i}(x);
                else
                    %...or with parameters...
                    x=execfun{i}(x,param{i});
                end
            else
                %...print message...
                %cprintf('TEXT',['\nWARNING in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\n'])
                fprintf(['\nWARNING in function ',funcName,':','\n',varName,' did not satisfy requirement ', dummyf.function,'\n'])
            end
        end
    end
      
    i=i+1;
    CONTINUE = ( i<=Ntests ) && ( RESULT(i-1)>=0 );
end



        
