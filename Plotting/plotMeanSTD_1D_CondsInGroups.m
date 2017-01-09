function [h, hsub] = plotMeanSTD_1D_CondsInGroups(mean,std,Measure, Colors, Groups, Conds, Ns,cfgMeasure,LineStyle)

if nargin<9
    LineStyle = '-';
end
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
    

h = figure('Name',titleName);
maxVal = -inf;
minVal = inf;
for iG=1:Ng;
    
    hsub(iG) = subplot(1,Ng,iG);
    hold on;
    
    for iC=1:Nc;
        
        iX = (iG-1)*Nc + iC;
        
        N=length(mean{iG,iC});
        
        x = [1:N].';
        
        ll = mean{iG,iC} - std{iG,iC};
        ul = mean{iG,iC} + std{iG,iC};
        
        minVal = min(minVal,min(ll));
        maxVal = max(maxVal,max(ul));
        
        p=patch([x; x(end:-1:1)],[ul; ll(end:-1:1)],ones(1,2*N),'edgecolor',Colors{iX},'facecolor',Colors{iX},'facealpha',0.25,'edgealpha',0);
        hasbehavior(p,'legend',false);
        plot(x,mean{iG,iC},'color',Colors{iX},'LineStyle',LineStyle,'linewidth',2);
    end
end

for iG=1:Ng;
    subplot(hsub(iG));
    hold on;
    title(Groups{iG})
    set(gca,'ylim',[minVal,maxVal]);
    legend(Conds)
    grid on;
    hold off
end

        