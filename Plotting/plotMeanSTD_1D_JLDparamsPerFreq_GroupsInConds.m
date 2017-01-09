function [h, hsub] = plotMeanSTD_1D_JLDparamsPerFreq_GroupsInConds(mean,std,Measure, Colors, Groups, Conds, Ns,cfgMeasure, Nm,iM)

freqs = cfgMeasure.freqs;

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

LineStyle = '-';%'none';
Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

%Make legend
for iC=1:Nc;
    xlbl{iC} = [Conds{iC}];
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

maxVal = -inf;
minVal = inf;

for iC=1:Nc;
    subplot(Nm,Nc,Nc*(iM-1)+iC);
    hold on;
    for iG=1:Ng;
        
        ll{iG,iC} = mean{iG,iC} - std{iG,iC};
        ul{iG,iC} = mean{iG,iC} + std{iG,iC};
        minVal = min(minVal,min(ll{iG,iC}));
        maxVal = max(maxVal,max(ul{iG,iC}));
        
        errorbar(x+0.5*(iG-3)-0.5,mean{iG,iC},std{iG,iC},'color',Colors{Nc*(iG-1)+iC},'LineStyle',LineStyle,'linewidth',2,'Marker','x','MarkerSize',5)
        
    end
end

for iC=1:Nc;
    subplot(Nm,Nc,Nc*(iM-1)+iC);
    hold on;
    title(xlbl{iC})
    set(gca,'xlim',[x(1)-1,x(end)+1]);
    set(gca,'ylim',[minVal,maxVal])
    grid on;
    set(gca,'xtick',x)
    if iM==Nm
        if iC==2
            xlabel('Frequency (Hz)')
        end
    else
        set(gca,'xticklabel',{})
    end
    if iC==1
        ylabel(Measure)
    else
        set(gca,'yticklabel',{})
    end
    if iM==1
        legend(Groups)
    end
    hold off
end