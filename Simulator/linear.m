function [dxdt p] = linear(x,p,t)

dxdt = (-x+p.x0)./p.tau; 