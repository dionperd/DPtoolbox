function h = plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin)

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

%Make legend
iX=0;
if (Ng==1)
    for iC=1:Nc;
        iX=iX+1;
        xlbl{iX} = [Conds{iC}];
    end
    Ndof = Nc;
elseif (Nc==1)
    for iG=1:Ng;
        iX=iX+1;
        xlbl{iX} = [Groups{iG}];
    end
    Ndof = Ng;
else
    for iG=1:Ng;
        for iC=1:Nc;
            iX=iX+1;
            xlbl{iX} = [Groups{iG},'-',Conds{iC}];
        end
    end
    Ndof = Ng*Nc;
end

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
title(titleName);
mean = mean';
std = std';
for iX=1:Ndof;
    hold on;
    bar(iX,mean{iX},0.5,'FaceColor',Colors{iX});
    hold on;
    errorbar(iX,mean{iX},mean{iX}-std{iX},mean{iX}+std{iX},'Color','k')
end
axis tight;
set(gca,'xlim',[0 Ndof+1],'xtick',[1:Ndof],'xticklabel',xlbl);
hold off
        