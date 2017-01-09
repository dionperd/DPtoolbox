function [R, G] = DPcalcAutoCorr(data,scales,num_tpts,num_vars,num_sc, norm, slopeInds) 

m = mean(data);
v = var(data);
maxlag = scales(end);
R = zeros(num_sc,num_vars);
for  iV = 1:num_vars;
    temp = xcorr( data(:,iV)-m(iV), maxlag) /v(iV);
    temp = temp(maxlag+2:end); %take scales after >0
    R(:,iV) = temp(scales);
    clear temp;
end
if norm
    R = R/num_tpts;
end

if ~isempty(slopeInds)
    for  iV = 1:num_vars;
        G(iV) = DPcalcSlope(log(scales),log(R(:,iV)),slopeInds);
    end
else
   G=zeros(1,num_vars);
end