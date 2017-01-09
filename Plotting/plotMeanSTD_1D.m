function [h, hsub] = plotMeanSTD_1D(mean,std,Measure, Colors, Groups, Conds, Ns,LineStyle,cfgMeasure)

if nargin<8
    LineStyle = '-';
end

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
mean = mean';
std = std';
maxVal = -inf;
minVal = inf;
xx = 1;
if isempty(cfgMeasure)
    xx=0;
else
    if ~isstruct(cfgMeasure)
        xx=0;
    else
        if ~isfield(cfgMeasure,'x')
            xx=0;
        end
    end
end
for iX=1:Ndof;
    
    N=length(mean{iX});
    mean{iX} = mean{iX}(:);
    std{iX} = std{iX}(:);
    ll{iX} = mean{iX} - std{iX}; 
    ul{iX} = mean{iX} + std{iX};
    minVal = min(minVal,min(ll{iX}));
    maxVal = max(maxVal,max(ul{iX}));
    
    hsub(iX) = subplot(Ng,Nc,iX);
    hold on;
    if xx==1
        x = cfgMeasure.x;
    else
        x = reshape(1:length(mean{iX}),size(mean{iX}));
    end
    p=patch([x; x(end:-1:1)],[ul{iX}; ll{iX}(end:-1:1)],ones(1,2*N),'edgecolor',Colors{iX},'facecolor',Colors{iX},'facealpha',0.25,'edgealpha',0);
    hasbehavior(p,'legend',false);
    plot(x,mean{iX},'color',Colors{iX},'LineStyle',LineStyle,'linewidth',2);
end

for iX=1:Ndof;
    subplot(hsub(iX));
    hold on;
    title(xlbl{iX})
    set(gca,'ylim',[minVal,maxVal]);
    grid on;
    hold off
end

        