function [h, hsub] = plotMeanSTD_2D_JLDparamsAlltogether(mean,std,Measure, Colors, Groups, Conds, Ns,cfgMeasure)

scales = cfgMeasure.scales;
freqs = cfgMeasure.freqs;

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

Nm = numel(Measure);

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
    titleName = [Measure,', subjects'' mean'];
    x=mean.';
else
    titleName = [Measure,', subjects'' mean / standard error'];
    %Calculate mean / standard error ratio
    for iG=1:Ng;
        for iC=1:Nc;
            x{iC,iG} = mean{iG,iC} ./ (std{iG,iC}/sqrt(Ns(iG))) ; 
        end
    end
end
    
for iM = 1:Nm;
    
    if isempty(Ns)
        titleName = [Measure{iM},', subjects'' mean'];
    else
        titleName = [Measure{iM},', subjects'' mean / standard error'];
    end
    
    h(iM) = figure('Name',titleName);%, 'Visible','off'
    
    maxVal = -inf;
    minVal = inf;
    for iX=1:Ndof;
        minVal = min(minVal,min(min(x{iX}(:,:,iM))));
        maxVal = max(maxVal,max(max(x{iX}(:,:,iM))));
    end
    
    for iX=1:Ndof;
        subplot(Ng,Nc,iX);
        hold on;
        title(xlbl{iX})
        imagesc(scales,freqs,squeeze(x{iX}(:,:,iM)).',[minVal,maxVal])
        set(gca,'Ydir','normal');
        xlabel('Time scales (sec)')
        ylabel('Frequency (Hz)')
        axis tight
        if iX==1:Ndof;
            colorbar;
        end
        hold off
    end
end
