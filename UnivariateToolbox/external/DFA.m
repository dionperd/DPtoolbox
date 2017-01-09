function output1 = DFA(DATA,win_length,order)

% DFA computes the Detrended Fluctuations of a given time series BB as a
% function of scales and delivers it as output1 of this function.
%
% DATA = timeseries
% win_length = desired window length of a given box, needs to be specified in fucntion  
% order = order of polynomial fit, for DFA to be set as 1
%
% 10/04/2011 - Original version: Viktor Jirsa
%              this version used the functions LinearFit and LinearLSQ
% 03/05/2011 - VJ: simplified version of the algorithm using the MATLAB
%              functions polyfit and polyval (found on the internet, tested
%              by VJ)


N=length(DATA);        % total length of time series
n=floor(N/win_length); % number of boxes on a given scale
N1=n*win_length; 
y=zeros(N1,1); Yn=zeros(N1,1); fitcoef=zeros(n,order+1);

mean1=mean(DATA(1:N1));
for i=1:N1
    
    y(i)=sum(DATA(1:i)-mean1);
end
y=y';
for j=1:n
    fitcoef(j,:)=polyfit(1:win_length,y(((j-1)*win_length+1):j*win_length),order);
end

for j=1:n
    Yn(((j-1)*win_length+1):j*win_length)=polyval(fitcoef(j,:),1:win_length);
end

output1=sqrt(sum((y'-Yn).^2)/N1);

return