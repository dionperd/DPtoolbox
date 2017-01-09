function [dxdt pout] = linearND(t,x,p)

pout=p;

dxdt = -p.a.*x+p.b;