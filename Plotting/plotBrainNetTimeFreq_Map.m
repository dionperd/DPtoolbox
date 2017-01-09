function [h, hsub] = plotBrainNetTimeFreq_Map(h,C,type,clim,t,tSel,tlbl,f,fSel,flbl,figTitle,varargin)

tic1=tic;
Nt = numel(tSel);
Nf = numel(fSel);

FigPos = [1 1 15 28];
if (h==0)
    h = figure('color','w','Units','centimeters', 'PaperPositionMode','auto','Position',FigPos);
else
    set(h,'color','w','Units','centimeters', 'PaperPositionMode','auto','Position',FigPos);
end
offsetWidth = 1;%cm
offsetHeight = 1;%cm
offsetWidth = offsetWidth/FigPos(3);
offsetHeight = offsetHeight/FigPos(4);
sbWidth = (1-offsetWidth)/Nt; 
sbHeight = (1-2*offsetHeight)/Nf;
marginHeight = 0.1*sbHeight;
marginWidth = 0.1*sbWidth;

setClim = iscell(clim);

for iF = 1:Nf;
    
    if setClim
        thisClim = [clim{1}(iF), clim{2}(iF)];
    else
        thisClim = clim;
    end
    for iT = 1:Nt;
             
        iF
        iT
        tic
        
        if strcmpi(type,'within')
            %select the time-frequency window to be plotted for each subject
            for iS=1:2;
                net{iS} = C{iS}(t(tSel{iT}),f(fSel{iF}),:,:);
                %take the mean over the time
                net{iS} = mean(net{iS},1);
                %take the mean over the frequencies
                net{iS} = mean(net{iS},2);
                
                net{iS}=squeeze(net{iS});
            end
            net{iS} = net{iS}.';
            net = net{1} + net{2};
           
        else          
            %select the time-frequency window to be plotted
            net = C(t(tSel{iT}),f(fSel{iF}),:,:);
            %take the mean over the time
            net = mean(net,1);
            %take the mean over the frequencies
            net = mean(net,2);
            
            net=squeeze(net);
            
           if ~strcmpi(type,'between')
                net = net.';
           end
        end
         
        if (iF==1)
            axXlbl = tlbl{iT};
        else
            axXlbl = '';
        end
        
        if (iT==1)
            axYlbl = flbl{iF};
        else
            axYlbl = '';
        end
        
        if (iF==Nf) && (iT==figTitle{1})
            axTitle=figTitle{2};
        else
            axTitle='';
        end
        
        hsub(iF,iT) = subplot('Position',[(iT-1)*sbWidth+marginWidth + offsetWidth, (iF-1)*sbHeight+marginHeight + offsetHeight, sbWidth-2*marginWidth, sbHeight-2*marginHeight]);%subplot(Nf,Nt,iSub);

        %For colormap connectivities:
        Nch = size(net,1);
        net(net==0) = thisClim(1)-(thisClim(2)-thisClim(1))/256 ;%clim(1);
        imagesc(1:Nch,1:Nch,net),%,clim);
        cmap = colormap(jet(256));
        cmap(2:end+1,:) = cmap;
        cmap(1,:) = [1 1 1];
        cmap(end+1,:) = [1 1 1];
        colormap(cmap)
        %colorbar;
        set(gca,'Ydir','reverse','xticklabel',{},'yticklabel',{},'clim',thisClim);%);
        box off;
        axis off;
        
%         ax=pcolor(net);
%         colormap(jet);
%         %  shading flat
%         % shading interp
%         set(ax,'edgecolor','none')
%         axis off;
%         axis ij;
%         
        xlabel(axXlbl); ylabel(axYlbl); title(axTitle);
        set(findall(gca, 'type', 'text'), 'visible', 'on');
        
        clear net;
        toc
        
    end
    
end
toc(tic1);

