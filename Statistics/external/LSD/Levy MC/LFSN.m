function X=LFSN(alpha,beta,H,a,b,m,M,N)
	T=m*(2*M+N-1);
	gv=g(a,b,[1:2*m*M]'/m-M,H-1/alpha);
    X=ifft(fft(gv,T).*fft(stblrnd(alpha,beta,(1/m)^(1/alpha),0,T,1)));
	X=real(X([1:N]*m));