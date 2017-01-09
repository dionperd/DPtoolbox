function [h, hsub] = plotMeanSTD_BrainNet(mean,std,Measure, Colors, Groups, Conds, Ns,varargin)

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

cfg = varargin{1};

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
    for iG=1:Ng;
        for iC=1:Nc;
            x{iC,iG} = mean{iG,iC} .' ; 
            x{iC,iG} = DPunPackBrainNet(mean,cfg);
        end
    end
else
    titleName = [Measure,', subjects'' mean / standard error'];
    %Calculate mean / standard error ratio
    for iG=1:Ng;
        for iC=1:Nc;
            x{iC,iG} = mean{iG,iC} ./ (std{iG,iC}/sqrt(Ns(iG))) ; 
            x{iC,iG} = DPunPackBrainNet(mean,cfg);
        end
    end
end


Nt = size(x{1,1},1);

for iT = 1:Nt;
    
for iF = 1:cfg.Nf;
    
    net = x{iC,iG}(
    h(iF) = figure('Name',[titleName,': f = ',num2str(cfg.fc(iF))]);
    
    for iX=1:Ndof;
        hsub(iF,iX) = subplot(Ng,Nc,iX);
        hold on;
        figure('color','w');
        plotBrainNet(C,hsub(iF,iX))
%         DPBrainNet_MapCfg();
%         hBN = gcf;
%         axBN = findall(hBN,'type','axes');
%         hsub(iF,iX)=copyobj(allchild(axBN(1)),hsub(iF,iX));
        colorbar;
        title(xlbl{iX})
        axis tight
        hold off
%         close(hBN);
%         delete hBN;
    end

end

end

title('test');
axis off;