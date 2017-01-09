function [ x, a, b ] = DPbutterFilt(x, fs, flim, type, order)

if nargin<5
    order = 4;
end

if nargin<4
    if length(flim)>1
        type = 'band';
    else
        type = 'low';
    end
end

%Normalize frequency:
fnyq = fs/2; %Nyquist frequency
Wn = flim/fnyq;

switch type
    
    case 'high'
        [b,a] = butter(ceil(order/2),Wn,'high');  
    case 'stop'
        [b,a] = butter(ceil(order/2),Wn,'stop');    
    otherwise
        [b,a] = butter(ceil(order/2),Wn);   
end
% freqz(b,a)
% dataIn = randn(1000,1);
x = filter(b,a,x);