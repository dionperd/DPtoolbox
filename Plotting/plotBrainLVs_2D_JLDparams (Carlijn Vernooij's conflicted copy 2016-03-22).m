function plotBrainLVs_2D_JLDparams(LV,BS,opts,Measure,scales,freqs)

if isempty(BS) 
    minVal = min(LV(:));
    maxVal = max(LV(:));
    imagesc(LV.',[minVal,maxVal])
else
    if strcmpi(opts.method,'BootsRatios')
       x=BS.'; 
    else
       x=LV.';
    end
    minVal = min(x(:));
    maxVal = max(x(:));
    alphadata = ones(size(BS));
    inds = abs(BS)<opts.th;
    minBS = min(BS(:));
    alphadata(inds) =  opts.Range(1) + (BS(inds) - minBS)*( diff(opts.Range)/(opts.th-minBS) );
    imagesc(scales,freqs,x,'AlphaData',alphadata.');
    set(gca,'clim',[minVal,maxVal]);
    xlabel('Time scales (sec)')
    ylabel('Frequency (Hz)')
end
%set(gca,'Ydir','normal');
axis tight
colorbar;
hold off

        