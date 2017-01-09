function [h, hsub] = plotMeanSTD_2D_JLDparams(mean,std,Measure, Colors, Groups, Conds, Ns,scales,freqs)

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
    

h = figure('Name',titleName);
maxVal = -inf;
minVal = inf;
for iX=1:Ndof;
    minVal = min(minVal,min(min(x{iX})));
    maxVal = max(maxVal,max(max(x{iX})));    
end

for iX=1:Ndof;
    subplot(Ng,Nc,iX);
    hold on;
    title(xlbl{iX})
    imagesc(scales,freqs,x{iX}.',[minVal,maxVal])
    set(gca,'Ydir','reverse');
    xlabel('Time scales (sec)')
    ylabel('Frequency (Hz)')
    axis tight
    %colorbar;
    hold off
end

        