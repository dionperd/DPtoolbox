function x=DProw2column(x)

%This function converts row wise matrices or vectors to column wise ones

%Input:
%x: 1D or 2D array of any class

%Output
%x: equals x' of the input if xSIZE(1)<xSIZE(2), where xSIZE=size(x) at the
%   input


%Create the input parser
p=inputParser;
p.FunctionName = 'DProw2column';
p.CaseSensitive=false; %NOT case sensitive
p.KeepUnmatched = false; %do not accept inputs undeclared here
p.StructExpand = false; %accept structures as single inputs

xSIZE = size(x);

p.addRequired('x',@(x)(length(xSIZE)==2));

p.parse(x);

if xSIZE(1)<xSIZE(2)
    x=x';
end




