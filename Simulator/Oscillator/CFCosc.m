function [dxdt pout] = CFCosc(t,x,p)


a(1,1) = p.a{1,1}(t);
a(1,2) = p.a{1,2}(t);
a(2,1) = p.a{2,1}(t);
a(2,2) = p.a{2,2}(t);
b(1) = p.b{1}(t);
b(2) = p.b{2}(t);
d(1) = p.d{1}(t);
d(2) = p.d{2}(t);


dxdt(1,1) = x(1) * ( 1 -a(1,1)*x(1)^2 -a(1,2)*x(3)^2 - b(1)*sin(p.n*x(2) - p.m*x(4)) );

dxdt(1,2) = p.omega(1) - d(1)*sin(p.n*x(2) - p.m*x(4)) ;

dxdt(1,3) = x(3) * ( 1 -a(2,1)*x(1)^2 -a(2,2)*x(3)^2 - b(2)*sin(p.m*x(4) - p.n*x(2)) );

dxdt(1,4) = p.omega(2) - d(2)*sin(p.m*x(4) - p.n*x(2)) ;



pout.omega = p.omega;
pout.m = p.m;
pout.n = p.n;
pout.a(1,1) =a(1,1);
pout.a(1,2) =a(1,2);
pout.a(2,1) =a(2,1);
pout.a(2,2) =a(2,2);
pout.b(1) =b(1);
pout.b(2) =b(2);
pout.d(1) =d(1);
pout.d(2) =d(2);


