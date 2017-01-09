function plotBrainLVs_1D(LV,BS,opts,Measure,LineColor,LineStyle,varargin)

if nargin<6
    LineStyle = '-';
end
if nargin<5
    LineColor = 'b';
end
if isempty(BS) 
    plot(LV,'color',LineColor,'LineStyle',LineStyle,'LineWidth',2);
else
    if strcmpi(opts.method,'BootsRatios')
        plot(BS,'color',LineColor,'LineStyle',LineStyle,'LineWidth',2);
        hold on;
        N = length(BS);
        LineTH = opts.th*ones(1,N);
        plot(LineTH,'color',LineColor,'LineStyle','--','LineWidth',1);
        plot(-LineTH,'color',LineColor,'LineStyle','--','LineWidth',1);

    else
        plot(LV,'color',LineColor,'LineStyle',LineStyle,'LineWidth',2);
        hold on;
        plot(LV(abs(BS)>=opt.th),'*','MarkerFaceColor',LineColor,'MarkerEdgeColor',LineColor,'MarkerSize',2);
    end

end
grid on;
ylabel(Measure)
axis tight


        