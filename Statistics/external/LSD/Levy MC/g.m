function g=g(a,b,s,d)
	A=gpower(s,d,1)-gpower(s-1,d,1);
	B=gpower(s,d,0)-gpower(s-1,d,0);
	g=a*A+b*B;