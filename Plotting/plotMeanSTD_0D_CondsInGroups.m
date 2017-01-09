function h = plotMeanSTD_0D_CondsInGroups(mean,std,Measure, Colors, Groups, Conds, Ns,varargin)

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

if isempty(Ns)
    titleName = [Measure,', subjects'' mean and standard deviation'];
else
    titleName = [Measure,', subjects'' mean and standard errors'];
    %Calculate standard error
    for iG=1:Ng;
        for iC=1:Nc;
            std0(iG,iC) = std{iG,iC}/sqrt(Ns(iG)); 
        end
    end
end

h = figure('Name',titleName);
title(titleName);
iX=0;
for iG=1:Ng;
    for iC=1:Nc;
        iX=iX+1;
        iClr = (iG-1)*Nc + iC;
        hold on;
        bar(iX,mean{iG,iC},0.5,'FaceColor',Colors{iClr});
        hold on;
        errorbar(iX,mean{iG,iC},mean{iG,iC}-std{iG,iC},mean{iG,iC}+std{iG,iC},'Color','k')
        xlbl{iX} = [Groups{iG},'-',Conds{iC}];
    end
    iX=iX+1;
    xlbl{iX} = '';
end
axis tight;
set(gca,'xlim',[0 iX],'xtick',[1:iX],'xticklabel',xlbl);
hold off
        