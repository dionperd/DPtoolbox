function x = DPsurrTrialShuf(x,D,Ntr)

%This function creates a surrogate set of data series by shuffling randomly
%the indexes of all columns except for the first one

%Inputs:
%   -x: the data series columnwise matrix
%   -D: the number of dimensions of the set
%   -Ntr: the number of trials


%Output:
%   -x: the surrogate data series


%For each dimension starting from the second one...
for iD=2:D;
    
    %...create a random set of indices...
    randIND = randperm(Ntr);
    
    %...and shuffle the data.   
    x(:,iD,:) = x(:,iD,randIND);
   
end
    