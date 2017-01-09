function [h, hsub] = plotMeanSTD_1D_JLDparamsPerFreq_CondsInGroups(mean,std,Measure, Colors, Groups, Conds, Ns,cfgMeasure, Nm,iM)

freqs = cfgMeasure.freqs;

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

LineStyle = '-';%'none';
Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

%Make legend
for iG=1:Ng;
    xlbl{iG} = [Groups{iG}];
end
Ndof = Ng*Nc;


if isempty(Ns)
    %titleName = [Measure,', subjects'' mean and standard deviation'];
else
    %titleName = [Measure,', subjects'' mean and standard errors'];
    %Calculate standard error
    for iG=1:Ng;
        for iC=1:Nc;
            std{iG,iC} = std{iG,iC}/sqrt(Ns(iG)); 
        end
    end
end
    
x=freqs;
mean = mean';
std = std';
maxVal = -inf;
minVal = inf;
iX=0;
for iG=1:Ng;
    subplot(Nm,Ng,Ng*(iM-1)+iG);
    hold on;
    for iC=1:Nc;
        iX = iX+1;
        N=length(mean{iX});
        
        ll{iX} = mean{iX} - std{iX};
        ul{iX} = mean{iX} + std{iX};
        minVal = min(minVal,min(ll{iX}));
        maxVal = max(maxVal,max(ul{iX}));
        
        
        %     p=patch([x; x(end:-1:1)],[ul{iX}; ll{iX}(end:-1:1)],ones(1,2*N),'edgecolor',Colors{iX},'facecolor',Colors{iX},'facealpha',0.25,'edgealpha',0);
        %     hasbehavior(p,'legend',false);
        %     plot(x,mean{iX},'color',Colors{iX},'LineStyle',LineStyle,'linewidth',2);
        
        errorbar(x+0.5*(iC-2),mean{iX},std{iX},'color',Colors{iX},'LineStyle',LineStyle,'linewidth',2,'Marker','x','MarkerSize',5)
        
    end
end

for iG=1:Ng;
    subplot(Nm,Ng,Ng*(iM-1)+iG);
    hold on;
    title(xlbl{iG})
    set(gca,'xlim',[x(1)-1,x(end)+1]);
    set(gca,'ylim',[minVal,maxVal])
    grid on;
    set(gca,'xtick',x)
    if iM==Nm
        if iG==2
            xlabel('Frequency (Hz)')
        end
    else
        set(gca,'xticklabel',{})
    end
    if iG==1
        ylabel(Measure)
    else
        set(gca,'yticklabel',{})
    end
    if iM==1
        legend(Conds)
    end
    hold off
end