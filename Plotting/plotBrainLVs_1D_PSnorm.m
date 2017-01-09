function plotBrainLVs_1D_PSnorm(LV,BS,opts,Measure,cfgMeasure)

if nargin<7
    MarkerStyle = '*';
end
if nargin<6
    Color = 'k';
end

freqs = cfgMeasure.freqs;
Nf = length(freqs);

if isempty(BS) 
    minVal = min(LV(:));
    maxVal = max(LV(:));
    maxAbsVal = max(minVal, maxVal);
    hold on;
    for iF = 1:Nf;
        plot(freqs(iF),LV(:,iF),'Marker',MarkerStyle,'markerfacecolor',Color,'MarkerEdgeColor',Color,'MarkerSize',2,'Linestyle','none');
    end
    xlabel('Frequency (Hz)')
    axis([freqs(1)-0.5, freqs(Nf)+0.5, minVal-0.1*maxAbsVal, maxVal+0.1*maxAbsVal])
    set(gca,'xtick',freqs)
else
    if strcmpi(opts.method,'BootsRatios')
       x=BS; 
    else
       x=LV;
    end
    minVal = min(x(:));
    maxVal = max(x(:));
    maxAbsVal = max(minVal, maxVal);
    hold on;
    if strcmpi(opts.method,'BootsRatios')
        for iF = 1:Nf;
            plot(freqs(iF),x(:,iF),'Marker',MarkerStyle,'markerfacecolor',Color,'MarkerEdgeColor',Color,'MarkerSize',2,'Linestyle','none');
        end
        plot([freqs(1)-0.5,freqs(Nf)+0.5],[opts.th, opts.th],'r','LineStyle','--','LineWidth',2)
        plot([freqs(1)-0.5,freqs(Nf)+0.5],[-opts.th, -opts.th],'r','LineStyle','--','LineWidth',2)
    else
        for iF = 1:Nf;
            plot(freqs(iF),x(abs(BS(:,iF))<opts.th,iF),'Marker',MarkerStyle,'markerfacecolor',Color,'MarkerEdgeColor',Color,'MarkerSize',2,'Linestyle','none');
        end
        for iF = 1:Nf;
            plot(freqs(iF),x(abs(BS(:,iF))>=opts.th,iF),'Marker',MarkerStyle,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',2,'Linestyle','none');
        end
    end
    axis([freqs(1)-0.5, freqs(Nf)+0.5, minVal-0.1*maxAbsVal, maxVal+0.1*maxAbsVal])
    set(gca,'xtick',freqs)
    xlabel('Frequency (Hz)')
end
hold off

        