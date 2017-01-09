function [h, hsub] = plotMeanSTD_1D_GroupsInConds_JLDparams(mean,std,Measure, Colors, Groups, Conds, Ns,scales)

LineStyle = '-';

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);


if isempty(Ns)
    titleName = [Measure,', subjects'' mean and standard deviation'];
else
    titleName = [Measure,', subjects'' mean and standard errors'];
    %Calculate standard error
    for iG=1:Ng;
        for iC=1:Nc;
            std{iG,iC} = std{iG,iC}/sqrt(Ns(iG)); 
        end
    end
end
    
x=scales;
h = figure('Name',titleName);
maxVal = -inf;
minVal = inf;
for iC=1:Nc;
    
    hsub(iC) = subplot(1,Nc,iC);
    hold on;
    
    for iG=1:Ng;
        
        iClr = (iG-1)*Nc + iC;
        
        N=length(mean{iG,iC});
        
        ll = mean{iG,iC} - std{iG,iC};
        ul = mean{iG,iC} + std{iG,iC};
        minVal = min(minVal,min(ll));
        maxVal = max(maxVal,max(ul));
        
        try
        p=patch([x; x(end:-1:1)],[ul; ll(end:-1:1)],ones(1,2*N),'edgecolor',Colors{iClr},'facecolor',Colors{iClr},'facealpha',0.25,'edgealpha',0);
        catch
            keyboard;
        end
        hasbehavior(p,'legend',false);
        plot(x,mean{iG,iC},'color',Colors{iClr},'LineStyle',LineStyle,'linewidth',2);
        
    end
end

for iC=1:Nc;
    subplot(hsub(iC));
    hold on;
    title(Conds{iC})
    set(gca,'ylim',[minVal,maxVal]);
    xlabel('Time scales (sec)')
    legend(Groups)
    grid on;
    hold off
end

        