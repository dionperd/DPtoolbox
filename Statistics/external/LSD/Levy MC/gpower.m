function g=gpower(s,d,type)
	g=s;
if (type==1) % positive
	g(s<0)=0;
	g(s>0)=s(s>0).^d;
else % negative
	g(s<0)=(-s(s<0)).^d;
	g(s>0)=0;
end