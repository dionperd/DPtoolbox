function x=DPsurrShufflTime(x,D,Ntr)

%This function shuffles randomly the points of the columns of x in order
%to create surrogate time series

%Get the initial size of x:
sizX0 = size(x);
reshapeFlag = length(sizX0)>2;
N=sizX0(1);

if reshapeFlag
    %Reshape so that only the first (time) dimension is retained and all
    %others are stacked to the second one
    x=reshape(x,[N,prod(sizX0(2:end))]);
end

x=x(randperm(N),:);
    
if reshapeFlag
    %Reshape to the original size:
    x=reshape(x,sizX0);
end
