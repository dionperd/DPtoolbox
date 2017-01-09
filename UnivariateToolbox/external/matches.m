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

n=length(y);
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
