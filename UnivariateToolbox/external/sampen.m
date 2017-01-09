function [e,C,M,R]=sampen(y,m,r,sflag,cflag)
%e=sampen(y,m,r)
%
%Input Parameters
%
%y  input signal vector
%m  maximum number of matches (default m=5)
%r  matching threshold (default r=.2)
%
%Output Parameters
%
%e  sample entropy calculations (m-1 values)
%
%Full usage:
%
%[e,C,M,R]=sampen(y,m,r,sflag,cflag)
%
%Input Parameters
%
%sflag    flag to standardize signal(default yes/sflag=1) 
%cflag    flag to use fast C code (default yes/cflag=1) 
%
%Output Parameters
%
%C        average number of m-matches per sample (m values)
%M        number of matches
%R        number of runs

if ~exist('m')|isempty(m),m=5;end
if ~exist('r')|isempty(r),r=.2;end
if ~exist('sflag')|isempty(sflag),sflag=1;end
if ~exist('cflag')|isempty(cflag),cflag=1;end
y=y(:);
n=length(y);
if sflag>0
   y=y-mean(y);
   s=sqrt(mean(y.^2));   
   y=y/s;
end
if cflag>0
   [M,R]=cmatches(y,n,r);
   M=double(M);
else
   [M,R]=matches(y,r);
end
k=length(M);
if k<m
   M((k+1):m)=0;
end
C=M(1:m)./(n-(1:m)');
%e=-diff(log(C));
e=zeros(m-1,1);
for i=1:m-1
   if C(i+1)>0
      e(i)=-log(C(i+1)/C(i));
   else
      e(i)=Inf;
   end
end

