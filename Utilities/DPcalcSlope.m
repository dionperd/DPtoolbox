function S = DPcalcSlope(x,y)

%This function finds the slope at every point i of the x-y curve, by fitting 
%a line 

%Inputs:
%-x: x coordinate, a vector of real numbers
%-y: y coordinate, a vector of real numbers
%-N: the length of x and y, positive integer

%Output:
%-S: the slope 


%...fit a line for those points...
p=polyfit(x, y, 1);

%...get the slope of the line
S=p(1);
