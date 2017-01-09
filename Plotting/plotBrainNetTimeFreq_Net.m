function [h, hsub] = plotBrainNetTimeFreq_Net(h,C,type,Measure,clim,t,tSel,tlbl,f,fSel,flbl,figTitle,varargin)

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
            Nch = size(net{1},1);
            Nch2 = 2*Nch;
            tempNet = zeros(Nch2,Nch2);
            tempNet(1:Nch,1:Nch) = net{1};
            tempNet((Nch+1):Nch2,(Nch+1):Nch2) = net{2};
            net = tempNet;
            clear tempNet;
        elseif strcmpi(type,'between')
            if strcmpi(Measure,'ICI')
                %select the time-frequency window to be plotted for each subject
                for iS=1:2;
                    net{iS} = C{iS}(t(tSel{iT}),f(fSel{iF}),:,:);
                    %take the mean over the time
                    net{iS} = mean(net{iS},1);
                    %take the mean over the frequencies
                    net{iS} = mean(net{iS},2);
                    
                    net{iS}=squeeze(net{iS});
                end
                Nch = size(net{1},1);
                Nch2 = 2*Nch;
                tempNet = zeros(Nch2,Nch2);
                tempNet(1:Nch,(Nch+1):Nch2) = net{1};
                tempNet((Nch+1):Nch2,1:Nch) = net{2};
                net = tempNet;
                clear tempNet;
            else
                %select the time-frequency window to be plotted
                net = C(t(tSel{iT}),f(fSel{iF}),:,:);
                %take the mean over the time
                net = mean(net,1);
                %take the mean over the frequencies
                net = mean(net,2);
                
                net=squeeze(net);
                
                Nch = size(net,1);
                Nch2 = 2*Nch;
                tempNet = zeros(Nch2,Nch2);
                tempNet(1:Nch,(Nch+1):Nch2) = net;
                net = tempNet;
                clear tempNet;
            end
        else
            %select the time-frequency window to be plotted
            net = C(t(tSel{iT}),f(fSel{iF}),:,:);
            %take the mean over the time
            net = mean(net,1);
            %take the mean over the frequencies
            net = mean(net,2);
            
            net=squeeze(net);
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
        
        %For BrainNet connectivities:
        ax = findall(hsub(iF,iT),'type','axes');
        %                                               surfFile,  nodeFile,  optionsFile
        plotBrainNet(net,ax,thisClim,axTitle,axXlbl,axYlbl,varargin{1},varargin{2},varargin{3});
        
        clear net;
        toc
        
    end
    
end
toc(tic1);

