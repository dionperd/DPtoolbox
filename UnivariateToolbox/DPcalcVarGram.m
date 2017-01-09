function [logV, VS] = DPcalcVarGram(x,scales,logScales,Nsc,slopeInds)

%This function calculates the Nsc scale variogram of time series x of N points
%according to V(sc) = mean( (x(i+sc) - x(i))^2 )
%VS: the slope of a logV-logScales diagram
%logV: the logarithm of the variogram


V=zeros(Nsc,1);
%For each scale...    
for iSc = 1:Nsc; %
    
    %...get it...
    sc=scales(iSc); 
    
    %...calculate V(sc)...
    V(iSc) = mean( ( x( (1+sc):end ) -  x( 1:end-sc ) ).^2 ); 
    
end

logV=log(V);

% VS=zeros(Nsc,1);
% if Nslope
%    inV = ~isinf(logV);
%     p=polyfit( logScales(inV), logV(inV), 1);
%     VS = p(1);
%    VS(inV) = DPcalcSlope(logScales(inV), logV(inV), sum(inV), Nslope);
% end

if ~isempty(slopeInds)
%     logScales=logScales(slopeInds);
%     logVtemp=logV(slopeInds);
%     inV = ~isinf(logVtemp) & ~isnan(logVtemp);
%     p=polyfit( logScales(inV), logVtemp(inV), 1);
%     VS=p(1);
    VS = DPcalcSlope(logScales,logV,slopeInds);
else
   VS=0;
end
logV=logV.';