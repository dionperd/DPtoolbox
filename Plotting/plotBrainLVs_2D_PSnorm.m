function plotBrainLVs_2D_PSnorm(LV,BS,opts,Measure,cfgMeasure)

freqs = cfgMeasure.freqs;

if isempty(BS) 
    minVal = min(LV(:));
    maxVal = max(LV(:));
    imagesc(LV.',[minVal,maxVal])
    Connections = size(LV,1);
else
    if strcmpi(opts.method,'BootsRatios')
       x=BS.'; 
    else
       x=LV.';
    end
    Connections = size(x,2);
    minVal = min(x(:));
    maxVal = max(x(:));
    alphadata = ones(size(BS));
    inds = abs(BS)<opts.th;
    minBS = min(BS(:));
    alphadata(inds) =  opts.Range(1) + (BS(inds) - minBS)*( diff(opts.Range)/(opts.th-minBS) );
    imagesc(Connections,freqs,x,'AlphaData',alphadata.');
    set(gca,'clim',[minVal,maxVal]);
    xlabel('Connections')
    ylabel('Frequency (Hz)')
end
%set(gca,'Ydir','normal');
axis tight
colorbar;
hold off

        