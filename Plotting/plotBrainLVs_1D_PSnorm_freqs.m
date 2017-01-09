function plotBrainLVs_1D_PSnorm_f(LV,BS,opts,Measure,cfgMeasure,Color,MarkerStyle)

if nargin<7
    MarkerStyle = '*';
end
if nargin<6
    Color = 'k';
end

f = cfgMeasure.f;
Nf = length(f);

if isempty(BS) 
    minVal = min(LV(:));
    maxVal = max(LV(:));
    maxAbsVal = max(minVal, maxVal);
    plot(f,LV,'markerfacecolor',Color,'markeredgecolor',Color,'Marker',MarkerStyle,'MarkerSize',5,'Linestyle','none');
    set(gca,'xtick',f)
    axis([f(1)-0.5, f(Nf)+0.5, minVal-0.1*maxAbsVal, maxVal+0.1*maxAbsVal])
else
    if strcmpi(opts.method,'BootsRatios')
        minVal = min(BS(:));
        maxVal = max(BS(:));
        maxAbsVal = max(minVal, maxVal);
        plot(f,BS,'markerfacecolor',Color,'MarkerEdgeColor',Color,'Marker',MarkerStyle,'MarkerSize',5,'Linestyle','none');
        hold on;
        LineTH = opts.th*ones(1,Nf);
        plot(f,LineTH,'color','r','LineStyle','--','LineWidth',2);
        plot(f,-LineTH,'color','r','LineStyle','--','LineWidth',2);
        axis([f(1)-0.5, f(Nf)+0.5, minVal-0.1*maxAbsVal, maxVal+0.1*maxAbsVal])
        set(gca,'xtick',f)
    else
        minVal = min(LV(:));
        maxVal = max(LV(:));
        maxAbsVal = max(minVal, maxVal);
        plot(f(abs(BS)<opts.th),LV(abs(BS)<opts.th),'Marker',MarkerStyle,'markerfacecolor',Color,'MarkerEdgeColor',Color,'MarkerSize',5,'Linestyle','none');
        hold on;
        plot(f(abs(BS)>=opts.th),LV(abs(BS)>=opts.th),'Marker',MarkerStyle,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'Linestyle','none');
        axis([f(1)-0.5, f(Nf)+0.5, minVal-0.1*maxAbsVal, maxVal+0.1*maxAbsVal])
        set(gca,'xtick',f)
    end

end
grid on;
ylabel(Measure)
xlabel('Frequency (Hz)')


        