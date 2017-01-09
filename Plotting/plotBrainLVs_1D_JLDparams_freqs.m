function plotBrainLVs_1D_JLDparams_freqs(LV,BS,opts,Measure,cfgMeasure,Color,MarkerStyle)

if nargin<7
    MarkerStyle = '*';
end
if nargin<6
    Color = 'b';
end

f = cfgMeasure.freqs;

if isempty(BS) 
    plot(f,LV,'markerfacecolor',Color,'markeredgecolor',Color,'Marker',MarkerStyle,'MarkerSize',5,'Linestyle','none');
else
    if strcmpi(opts.method,'BootsRatios')
        plot(f,BS,'markerfacecolor',Color,'MarkerEdgeColor',Color,'Marker',MarkerStyle,'MarkerSize',5,'Linestyle','none');
        hold on;
        N = length(BS);
        LineTH = opts.th*ones(1,N);
        plot(f,LineTH,'color',Color,'LineStyle','--','LineWidth',1);
        plot(f,-LineTH,'color',Color,'LineStyle','--','LineWidth',1);

    else
        plot(f(abs(BS)<opt.th),LV(abs(BS)<opt.th),'Marker',MarkerStyle,'markerfacecolor',Color,'MarkerEdgeColor',Color,'MarkerSize',5,'Linestyle','none');
        hold on;
        plot(f(abs(BS)>=opt.th),LV(abs(BS)>=opt.th),'Marker',MarkerStyle,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'Linestyle','none');
    end

end
grid on;
ylabel(Measure)
xlabel('Frequency (Hz)')
axis tight


        