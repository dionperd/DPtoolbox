function plotBrainLVs_1D_EMG_JLDparams(LV,BS,opts,Measure,scales,Color,LineStyle)

if nargin<7
    LineStyle = '-';
end
if nargin<6
    Color = 'b';
end
if isempty(BS) 
    plot(scales,LV,'color',Color,'LineStyle',LineStyle,'LineWidth',2);
else
    if strcmpi(opts.method,'BootsRatios')
        plot(scales,BS,'color',Color,'LineStyle',LineStyle,'LineWidth',2);
        hold on;
        N = length(BS);
        LineTH = opts.th*ones(1,N);
        plot(scales,LineTH,'color',Color,'LineStyle','--','LineWidth',1);
        plot(scales,-LineTH,'color',Color,'LineStyle','--','LineWidth',1);

    else
        plot(scales,LV,'color',Color,'LineStyle',LineStyle,'LineWidth',2);
        hold on;
        plot(scales,LV(abs(BS)>=opt.th),'*','MarkerFaceColor',Color,'MarkerEdgeColor',Color,'MarkerSize',2);
    end

end
grid on;
ylabel(Measure)
xlabel('Cycles')
axis tight


        