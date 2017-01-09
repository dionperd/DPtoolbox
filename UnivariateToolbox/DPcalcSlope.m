function S = DPcalcSlope(x,y,inds)

x=x(inds);
ytemp=y(inds);
inds = ~isinf(ytemp) & ~isnan(ytemp);
p=polyfit( x(inds), y(inds), 1);
S=p(1);