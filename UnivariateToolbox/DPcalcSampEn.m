function SE=DPcalcSampEn(y,m,r,n,sflag)

%This function calculates Sample Entropy using the Goldberger sampen code
%(see below also for inputs/outputs)

%--------------DP modified, see external folder for original core----------
%x=x(:);
%n=length(x);
if sflag>0
   y=zscore(y);
end


%Find matches:

% if cflag>0
%    [M,R]=cmatches(x,n,r);
%    M=double(M);
% else
%   [M,R]=matches(x,r);
%end
%k=length(M);

R=zeros(n-1,1);
for j=1:(n-1)
   i=1:(n-j);    
   d=abs(y(i+j)-y(i));
   a=d<r;
   a1=[0;a];
   a2=[a;0];
   aa=find(a1~=a2);
   kk=length(aa)/2;   
   rr=diff(reshape(aa,2,kk));
   for i=1:kk
      k=rr(i);
      R(k)=R(k)+1;
%       if k>50
%          [j k aa(2*i-1) aa(2*i)]
%       end      
   end      
end
k=find(R>0, 1, 'last' );
R=R(1:k);
M=zeros(k,1);
for i=1:k;
   M(i)=R(i);
   for j=(i+1):k;
      M(i)=M(i)+(j+1-i)*R(j);
   end
end   

%Calculate sample entropy
if k<m
   M((k+1):m)=0;
end
C=M(1:m)./(n-(1:m)');
%SE=-diff(log(C));
SE=zeros(m-1,1);
for i=1:m-1
   if C(i+1)>0
      SE(i)=-log(C(i+1)/C(i));
   else
      SE(i)=Inf;
   end
end


%Original functions not used here but copied-pasted in the code above:

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





function [M,R]=matches(y,r);
%[M,R]=matches(y,r)
%
%Input Parameters
%
%y  input signal vector
%r  matching threshold (default r=.2)
%
%Output Parameters
%
%M  number of matches
%R  number of runs

%n=length(y);
R=zeros(n-1,1);
for j=1:(n-1)
   i=1:(n-j);    
   d=abs(y(i+j)-y(i));
   a=d<r;
   a1=[0;a];
   a2=[a;0];
   aa=find(a1~=a2);
   kk=length(aa)/2;   
   rr=diff(reshape(aa,2,kk));
   for i=1:kk
      k=rr(i);
      R(k)=R(k)+1;
      if k>50
         [j k aa(2*i-1) aa(2*i)]
      end      
   end      
end
k=max(find(R>0));
R=R(1:k);
M=zeros(k,1);
for i=1:k
   M(i)=R(i);
   for j=(i+1):k
      M(i)=M(i)+(j+1-i)*R(j);
   end
end   
