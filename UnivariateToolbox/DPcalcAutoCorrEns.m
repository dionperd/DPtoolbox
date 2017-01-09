function [R, G] = DPcalcAutoCorrEns(data,scales,num_tpts,num_vars,num_sc, norm, slopeInds) 

m = mean(data);
v = var(data);
maxlag = scales(end);
R = zeros(num_sc,1);
for  iV = 1:num_vars;
    temp = xcorr( data(:,iV)-m(iV), maxlag) /v(iV);
    temp = temp(maxlag+2:end); %take scales after >0
    R(:,iV) = temp(scales);
    clear temp;
end
if norm
    R = R/num_tpts;
end
R = mean(R,2);

if ~isempty(slopeInds)
   G = DPcalcSlope(log(scales),log(R),slopeInds);
else
   G=0;
end