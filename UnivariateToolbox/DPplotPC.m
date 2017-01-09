function h=DPplotPC(x,phi,DphiCell,PCcell,statsResCell,cfg,method,measInds, Nmeasures, animation,filename)


%  Inputs:
%  -x: the time series of N time points, 
%      a matrix either Nx2 (for a single trial) or Nx2xNtr of Ntr trials data,
%      either real valued signal, or phase signal in  the interval [0 2pi]
%  -phi: phases of the signals wrapped into the interval [0 2pi)
%  -PCcell: results cell for the K different nm calculations 
%       of structures with fields with the names of the calculated measures
%       of size depending on the method:
%       1. vector 1 x Ntr for 'method='trial'
%       2. matrix Nout x Ntr for 'method='trialTime'
%       3. vector Nout x 1  for method='ensemble'    
%       4. vector Nout x 1 for 'method='ensembleTime', 
%          where Nout is the number of the window step calculations
% if method == 'trial' or 'trialTime', then mean values across trials are
% given in a substructure PC.trialMean as:
%       1. scalar for 'method='trial'
%       2. vector Nout(1) x 1 for 'trialWin'
%  -DphiCell: cell of the phase difference in the interval [-pi pi], a matrix of
%         dimensions Ncut x Ntr, for the K different nm calculations
%  -cfg: configuration structure corrected and complete
%  -method: 1,2,or 3 depending on whether cfg.method =  'trial', 
%           'trialTime','ensemble' respectively 
%  -measInds: the indexes of the measures to be calculated
%  -Nmeasures: how many measures are to be calculated
%  -thisfunCommands: cell of the command strings for the calculation of
%                    the selected measures


% Output:
% -h: Array of handles of the figures



%A usefull constant
TwoPi=2*pi;

%Cell of names of result's structure
MeasNames={'PCI','NCI','ACI','ICI12','ICI21','','';...
           'PLV','','','','','','';...
           'PPC','','','','','','';...
           'PLI','','','','','','';...
           'SE','','','','','','';...
           'CP','','','','','','';...
           'MI','','','','','','';...
           'PRPL','PRcn1','PRcn2','PRcn','PRcd1','PRcd2','PRcd';...
           };
NmeasPmeas = [5 1 1 1 1 1 1 7]; %Number of 'measures per measures' categories
NsubMeasures = sum(NmeasPmeas(measInds)); %Total number of measures calculated, submeasures included

%Measures' colors cell
MeasColors = {[0 0.4 0.2],[0.4 0.8 0.6],[0.4 0.6 0],[0 0.6 0],[0 1 0],[],[];...
              [0 0 1],[],[],[],[],[],[];...
              [0.2 0.4 0.8],[],[],[],[],[],[];...
              [0 1 1],[],[],[],[],[],[];...
              [0.6 0 1],[],[],[],[],[],[];...
              [0.8 0 0.8],[],[],[],[],[],[];...
              [0 0 0.6],[],[],[],[],[],[];...
              [0.4 0.2 0],[1 0.8 0],[0.8, 0.6 0.2],[0.6 0.4 0],[0.8 0.4 0],[1 0.4 0],[0.6 0.2 0];...
              };
SignalColors={'b','g'};
timeAxis = [cfg.time(1)-cfg.Ts cfg.time(end)+cfg.Ts];


histPhiTick = cfg.gridPhiEdgs; 
histDphiTick = histPhiTick - pi; 

histPhiTickLabelN=num2cell(histPhiTick/pi*180,[1,cfg.Nbins+1]);
histPhiTickLabelM=num2cell(histPhiTick/pi*180,[1,cfg.Nbins+1]);
histDphiTickLabel=num2cell(histDphiTick/pi*180,[1,cfg.Nbins+1]);
for ii=1:cfg.Nbins+1;
    histPhiTickLabelN{ii} = ['n',num2str(histPhiTickLabelN{ii}),'^{\circ}'];
    histPhiTickLabelM{ii} = ['m',num2str(histPhiTickLabelM{ii}),'^{\circ}'];
    histDphiTickLabel{ii} = [num2str(histDphiTickLabel{ii}),'^{\circ}'];   
end

if cfg.Ntr>1
    edgealpha = 0.1;
else
    edgealpha = 1;
end

scrsz = get(0,'ScreenSize');
fx=1;fy=scrsz(4)/2;fwidth=scrsz(3)/2.1;fheight=scrsz(4);



[K dummy]=size(cfg.nm); %number of different nm calculations

%Initialize hout
hout=nan(1,K);

for iNM=1:K;
    
    
    %Prepare data to be plotted:
    
    %This nm:
    nm=cfg.nm(iNM,:);
    TwoPiNM(1)=nm(1)*TwoPi;
    TwoPiNM(2)=nm(2)*TwoPi;
    
    %Get the caracters for plotting
    n=num2str(nm(1));
    m=num2str(nm(2));
    
    for iP=1:2;
        %Wrap phi in multiples of 2*pi
        phiNM(:,iP,:)=mod(nm(iP)*phi(:,iP,:),TwoPi);
    end
    
    %Get these Dphi and PC
    Dphi = DphiCell{nm(1),nm(2)};
    PC = PCcell{nm(1),nm(2)};
    
    
    %Create figure with subplots:
    h=figure('Position',[fx fy fwidth fheight],'color','w');
    %Subplots' structure
    sbp=[5,6];
    Nsubs=prod(sbp);
    %Where the measures are plotted
    sbpPos={1:sbp(2),sbp(2)+1:2*sbp(2),2*sbp(2)+1:3*sbp(2),Nsubs-2*sbp(2)+1:Nsubs};

            
    %Plot all trials' raw data and means commonly for all methods
    %Plot for each trial
    for iT=1:cfg.Ntr;
        
        for iP=1:2;            
            %Plot x
            hsubs(1)=subplot(sbp(1),sbp(2),sbpPos{1});
            hold on;
            hx(iP,iT)=patchline(cfg.time,x(:,iP,iT),'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',edgealpha,'linewidth',0.5);
            
            %Plot phi
            hsubs(2)=subplot(sbp(1),sbp(2),sbpPos{2});
            hold on;
            hphi(iP,iT)=patchline(cfg.timeCut,phiNM(:,iP,iT)/TwoPi,'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',edgealpha,'linewidth',0.5);
             
        end
        
        %Plot Dphi
        hsubs(3)=subplot(sbp(1),sbp(2),sbpPos{3});
        hold on;
        hDphi(iT)=patchline(cfg.timeCut,Dphi(:,iT)/pi,'edgecolor','k','edgealpha',edgealpha,'linewidth',0.5);
    end
    
    
    %Plot mean trajectories and label plots:
    
    %x
    xmax = max(x(:));
    xmax= xmax+0.01*abs(xmax);
    xmin = min(x(:));
    xmin= xmin-0.01*abs(xmin);
    subplot(hsubs(1))
    hold on;
    if cfg.Ntr>1
        for iP=1:2;
            hxM(iP)=patchline(cfg.time,mean(x(:,iP,:),3),'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',1,'linewidth',2);
        end
    end
    xlim(timeAxis)
    ylim([xmin xmax])
    ylabel('Signals')
    title([n,':',m,' phase coupling'])
    %legend('x_1', 'x_2')
    legend('Signal 1', 'Signal 2')
    hold off;
    
    %phi
    subplot(hsubs(2))
    hold on;
    if cfg.Ntr>1
        for iP=1:2;
            hphiM(iP)=patchline(cfg.timeCut,mean(phiNM(:,iP,:),3)/TwoPi,'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',1,'linewidth',2);
        end
    end
    xlim(timeAxis)
    ylim([-0.01 1.01 ])
    title([n,':',m,' phases'])
    %legend('\phi_1','\phi_2');
    set(gca,'ytick',[0:1/8:1]);%,'yticklabel',{'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','2*pi'}
    DPformat_ticks_ext_plus(gca,[],{'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'});
    grid on;
    hold off;
    
    %Dphi
    subplot(hsubs(3))
    hold on;
    if cfg.Ntr>1
        meanDPhi=mean(Dphi,2);
        hDphiM=patchline(cfg.timeCut,meanDPhi/pi,'edgecolor','r','facecolor','r','edgealpha',1,'linewidth',2);
    end
    xlim(timeAxis)
    ylim([-1.01 1.01 ])
    xlabel('Time (sec)')
    title([n,'\phi_1-',m,'\phi_2'])
    %set(gca,'ytick',[0:1/8:1]);
    set(gca,'ytick',[-1:1/4:1]);
    DPformat_ticks_ext_plus(gca,[],{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
    grid on;
    hold off;
    
    
    
    switch method
        
        case 1 %trial
            
            %Plot all trials' measure results
            hsubs(4)=subplot(sbp(1),sbp(2),sbpPos{4});
            hold on
            iiM=0;
            %For each measure we have calculated...
            for iM = 1:Nmeasures;
                %...and for each of the sumbmeasures of this measure...
                for iSubM = 1:NmeasPmeas(measInds(iM));
                    %...plot
                    iiM=iiM+1;
                    plot(iiM,PC.(MeasNames{measInds(iM),iSubM}),'.','color',MeasColors{measInds(iM),iSubM}(1,:),'markersize',15);
                    hold on;
                end
            end
            
            %Plot measures' means per trial
            iiM=0;
            xticklabel={};
            if cfg.Ntr>1
                %For each measure we have calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...plot
                        iiM=iiM+1;
                        hold on;
                        plot(iiM,PC.trialMean.(MeasNames{measInds(iM),iSubM}),'*','color',MeasColors{measInds(iM),iSubM}(1,:),'markersize',10);
                        xticklabel=[xticklabel,MeasNames{measInds(iM),iSubM}];
                    end
                end
            else %Or just place the xticklabel
                %For each measure we have calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...plot
                        iiM=iiM+1;
                        plot(iiM,PC.(MeasNames{measInds(iM),iSubM}),'*','color',MeasColors{measInds(iM),iSubM}(1,:),'markersize',10);
                        xticklabel=[xticklabel,MeasNames{measInds(iM),iSubM}];
                    end
                end
            end
            xlim([0 iiM+1])
            ylim([-1.01 1.01]);%'auto'
            set(gca,'xtick',[1:iiM],'xticklabel',xticklabel,'ytick',[-1:0.1:1]);
            grid on;
            xlabelM = 'Phase coupling measures';
            hxlabel=xlabel(xlabelM);
            hold off;
            
            
        case 2 %trialTime
            %Plot trials'means measure results
            hsubs(4)=subplot(sbp(1),sbp(2),sbpPos{4});
            if cfg.Nout(1)<=100
                marker = '-*';
            else
                marker = '-';   
            end
            hold on
            iiM=0;
            measLegend={};
            if cfg.Ntr>1
                %For each measure we have calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...plot
                        iiM=iiM+1;
                        plot(cfg.timeOut,PC.trialMean.(MeasNames{measInds(iM),iSubM}),marker,'color',MeasColors{measInds(iM),iSubM}(1,:),'markersize',5);
                        hold on;
                        measLegend=[measLegend;MeasNames{measInds(iM),iSubM}];
                    end
                end
            else
                %For each measure we have calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...plot
                        iiM=iiM+1;
                        hold on;
                        plot(cfg.timeOut,PC.(MeasNames{measInds(iM),iSubM}),marker,'color',MeasColors{measInds(iM),iSubM}(1,:),'markersize',5);
                        measLegend=[measLegend;MeasNames{measInds(iM),iSubM}];
                    end
                end
            end
            xlim(timeAxis)
            ylim([-1.01 1.01]);%'auto'
            set(gca,'xtick',cfg.timeOut,'ytick',[-1:0.1:1]);
            legend(measLegend)
            ylabel('Phase coupling measures')
            xlabelM = 'Time (sec)';
            hxlabel=xlabel(xlabelM);
            grid on;
            hold off;
            
            
        otherwise %ensemble
            
            if (cfg.Nwin==1)
                
                %Plot all trials' measure results
                hsubs(4)=subplot(sbp(1),sbp(2),sbpPos{4});
                hold on
                iiM=0;
                measLegend={};
                %For each measure we have calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...plot
                        iiM=iiM+1;
                        plot(cfg.timeOut,PC.(MeasNames{measInds(iM),iSubM}),'color',MeasColors{measInds(iM),iSubM}(1,:));
                        hold on;
                        measLegend=[measLegend;MeasNames{measInds(iM),iSubM}];
                    end
                end
                xlim(timeAxis)
                ylim([-1.01 1.01]);%'auto'
                yMlim = get(gca,'ylim');
                set(gca,'ytick',[-1:0.1:1])
                legend(measLegend)
                ylabel('Phase coupling measures')
                xlabelM = 'Time (sec)';
                xlabel(xlabelM)
                grid on;
                hold off;
                
            else
                
                %Plot all trials' measure results
                hsubs(4)=subplot(sbp(1),sbp(2),sbpPos{4});
                if cfg.Nout(1)<=100
                    marker = '-*';
                else
                    marker = '-';
                end
                hold on
                iiM=0;
                measLegend={};
                %For each measure we have calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        %...plot
                        iiM=iiM+1;
                        plot(cfg.timeOut,PC.(MeasNames{measInds(iM),iSubM}),marker,'color',MeasColors{measInds(iM),iSubM}(1,:),'markersize',5);
                        hold on;
                        measLegend=[measLegend;MeasNames{measInds(iM),iSubM}];
                    end
                end
                xlim(timeAxis)
                ylim([-1.01 1.01]);%'auto'
                yMlim = get(gca,'ylim');
                set(gca,'xtick',cfg.timeOut,'ytick',[-1:0.1:1]);
                legend(measLegend)
                ylabel('Phase coupling measures')
                xlabelM = 'Time (sec)';
                xlabel(xlabelM)
                grid on;
                hold off;
            end
    end %plot all trials' measure results
    
    %Save figure without animation
    saveas(h,[filename,'_',n,':',m,'.fig']);
    saveas(h,[filename,'_',n,':',m,'.png']);
    
    
    if any(animation>0)
        
        %Plot every k frames for the animation
        k=animation;
        
        switch method
            case 1
                                
                %Trial range selection
                iTrange = 1:k:cfg.Ntr;
                

                %Fade out mean trajectories
                if cfg.Ntr>1
                    for iP=1:2;
                        set(hxM(iP),'edgealpha',0.25,'facealpha',0.25);
                        set(hphiM(iP),'edgealpha',0.25,'facealpha',0.25);
                    end
                    set(hDphiM,'edgealpha',0.25,'facealpha',0.25);
                end
                maxHistY=cfg.Ncut; %maximum y value for the histogram plot
                
            case 2
 
                %Trial range selection
                iTrange = 1:k:cfg.Ntr;
                

                %Fade out mean trajectories
                if cfg.Ntr>1
                    for iP=1:2;
                        set(hxM(iP),'edgealpha',0.25,'facealpha',0.25);
                        set(hphiM(iP),'edgealpha',0.25,'facealpha',0.25);
                    end
                    set(hDphiM,'edgealpha',0.25,'facealpha',0.25);
                end
                maxHistY=cfg.Nwin; %maximum y value for the histogram plot
                
                if cfg.Nout(1)<=100
                    marker = '-o';
                else
                    marker = '-';
                end
                
                
            otherwise
                maxHistY=cfg.Ntr*cfg.Nwin; %maximum y value for the histogram plot
                %Time point range selection
                iTrange = 1:k:cfg.Ncalc;
                if cfg.Nwin>1
                    %We are going to need these for the window square plotting
                    C=ones(1,4);
                    syx = [xmin xmin xmax xmax];
                    sy = [-1 -1 1 1];
                    syM = [yMlim(1) yMlim(1) yMlim(2) yMlim(2)];
                end
                
        end
        %Calculation of histogram's yticks
        maxHistYs = [20 40 60 80 100 200 300 400 600 800 1000 1200 1600,...
                     2000 3000 4000 6000 8000 10000 12000 16000 20000,...
                     24000 28000 32000 30000 36000 40000 48000 60000,...
                     80000 100000 120000 160000 200000 240000 280000,...
                     300000 320000 360000 400000];
        maxHistYs=maxHistYs(maxHistYs>=maxHistY);
        maxHistY=maxHistYs(1);
        HistY=[0:maxHistY/4:maxHistY];
        
        
        %Rearrange figure
        
        %New subplot structure
        if (method==2)
            if cfg.Ncalc<=5
                sbp=[8,cfg.Ncalc];
                iSrange=1:cfg.Ncalc;
            else
                sbp=[8,5];
                iSrange = round(cfg.Ncalc*[2:4]/5);
                iSrange=[1 iSrange cfg.Ncalc];
            end
        else
            sbp=[7,6];
        end
        Nsubs=prod(sbp);
        sbpPos = {1:sbp(2),sbp(2)+1:2*sbp(2),2*sbp(2)+1:3*sbp(2),Nsubs-2*sbp(2)+1:Nsubs};
        
        %Make a temporary figure
        tempFig=figure('Position',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)]);
        
        for iF=1:4;
            %For each subplot...
            figure(tempFig); % Make sure tempFig is active
            %...create a dummy subplot...
            tempSub = subplot(sbp(1), sbp(2), sbpPos{iF});
            %...and get its position...
            tempPos = get(tempSub,'pos');
            %return to main figure...
            figure(h)
            %...change the position
            set(hsubs(iF),'pos',tempPos);
        end
        delete(tempFig)
        %Return to the figure
        figure(h)
        
        
        %Plot selected instants for the animation
        
        %         % Get figure size
        %         pos = get(gcf, 'Position');
        %         width = pos(3); height = pos(4);%
        %         % Preallocate data (for storing frame data)
        %         mov = zeros(height, width, 1, cfg.Ntr, 'uint8');
        
        
        %For each selected instant...
        Nframes = length(iTrange);
        iF=0; %frame meter and index
        for iT=iTrange;
            
            %Increase frame meter and index
            iF=iF+1;
            
            switch method
                case 1 %trial
                    thisDphi = Dphi(:,iT);
                    for ii=1:2
                        thisPhiNM(:,ii) = phiNM(:,ii,iT);
                    end
                case 3 %ensemble
                    thisDphi = squeeze( Dphi(cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),:) );
                    thisDphi=thisDphi(:);
                    for ii=1:2
                        temp = squeeze( phiNM(cfg.timePointWinEdgs(iT,1):cfg.timePointWinEdgs(iT,2),ii,:) );
                        thisPhiNM(:,ii) =temp(:);
                    end
            end
            
            if (method==2)
                iSub=0;
                for iS=iSrange;
                    iSub=iSub+1;
                    
                    for ii=1:2
                        thisPhiNM(:,ii) = phiNM(cfg.timePointWinEdgs(iS,1):cfg.timePointWinEdgs(iS,2),ii,iT);  
                    end
                    thisDphi = Dphi(cfg.timePointWinEdgs(iS,1):cfg.timePointWinEdgs(iS,2),iT);
                    
                    %Make compass (arrow) plot of Dphis
                    subplot(sbp(1),sbp(2),3*sbp(2)+iSub)
                    hold off;
                    DPcompass_transp(cos(thisDphi),sin(thisDphi),0,0.1,'-',gca);
                    %xlabel([n,'\phi_1-',m,'\phi_2'])
                    hold off;
                    
%                     %Make rose plot (circular histogram) of Dphis
%                     subplot(sbp(1),sbp(2),4*sbp(2)+iSub)
%                     hold off
%                     DProse(thisDphi,cfg.gridDphi);
%                     %xlabel([n,'\phi_1-',m,'\phi_2'])
%                     hold off;
                    
                    %Make histogram of Dphis
                    subplot(sbp(1),sbp(2),4*sbp(2)+iSub)
                    hold off;
                    
                    hist(thisDphi,cfg.gridDphi,'k');
                    hh = findobj(gca,'Type','patch');
                    set(hh,'FaceColor','k','EdgeColor','k')
                    xlim([-1.01 1.01 ]*pi)
                    ylim([-0.01 maxHistY])
                    title([n,'\phi_1-',m,'\phi_2']);hold on;
                    set(gca,'xtick',histDphiTick,'ytick',HistY);
                    xlabel([n,'\phi_1-',m,'\phi_2'])
                    DPformat_ticks_ext_plus(gca,histDphiTickLabel);
                    grid on;
                    hold off;

                    %Make histogram of phiNM
                    subplot(sbp(1),sbp(2),5*sbp(2)+iSub)
                    hold off;
                    NpointsPcalc=size(thisPhiNM,1);
                    p = hist3(thisPhiNM,{cfg.gridPhi cfg.gridPhi})/NpointsPcalc; 
                    imagesc(cfg.gridPhi,cfg.gridPhi,p.',[0 2/cfg.Nbins]);
                    axis square
                    Cmap=colormap(gray(100));
                    Cmap(1:100,:) = Cmap(100:-1:1,:);
                    colormap(Cmap);
                    %colorbar
                    set(gca,'YDir','normal')
%                     hold on;
%                     plot(thisPhiNM(:,1),thisPhiNM(:,2),'k.','markersize',0.5);
                    xlim([-0.01 1.01 ]*TwoPi)
                    ylim([-0.01 1.01 ]*TwoPi)
                    set(gca,'xtick',histPhiTick,'ytick',histPhiTick);
                    xlabel(['\fontsize{8}\bf',num2str(cfg.timeWinEdgs(iS,1)),' - ',num2str(cfg.timeWinEdgs(iS,2)), ' sec'],'color','k');
                    DPformat_ticks_ext_plus(gca,histPhiTickLabelN, histPhiTickLabelM);
                    grid on;
                    hold off;
                    
                    clear thisPhiNM thisDphi;
                    
                    
                end
                
            else
                %Make compass (arrow) plot of Dphis
                subplot(sbp(1),sbp(2),[3*sbp(2)+1,3*sbp(2)+2,4*sbp(2)+1,4*sbp(2)+2])
                hold off;
                DPcompass_transp(cos(thisDphi),sin(thisDphi),0,0.1,'-',gca);
                xlabel([n,'\phi_1-',m,'\phi_2'])
                hold off;
                
%                 %Make rose plot (circular histogram) of Dphis
%                 subplot(sbp(1),sbp(2),[3*sbp(2)+3,3*sbp(2)+4,4*sbp(2)+3,4*sbp(2)+4])
%                 hold off
%                 DProse(thisDphi,cfg.gridDphi);
%                 xlabel([n,'\phi_1-',m,'\phi_2'])
%                 hold off;
                
                %Make histogram of Dphis
                %subplot(sbp(1),sbp(2),[3*sbp(2)+5,3*sbp(2)+6,4*sbp(2)+5,4*sbp(2)+6])
                subplot(sbp(1),sbp(2),[3*sbp(2)+3,3*sbp(2)+4,4*sbp(2)+3,4*sbp(2)+4])
                hold off;
                hist(thisDphi,cfg.gridDphi,'k');
                hh = findobj(gca,'Type','patch');
                set(hh,'FaceColor','k','EdgeColor','k')
                xlim([-1.01 1.01 ]*pi)
                ylim([-0.01 maxHistY])
                %xlabel([n,'\phi_1-',m,'\phi_2'])
                set(gca,'xtick',histDphiTick,'ytick',HistY);
                xlabel([n,'\phi_1-',m,'\phi_2'])
                DPformat_ticks_ext_plus(gca,histDphiTickLabel);
                grid on;
                hold off;
                
                %Make histogram of phiNM
                %subplot(sbp(1),sbp(2),[3*sbp(2)+7,3*sbp(2)+8,4*sbp(2)+7,4*sbp(2)+8])
                subplot(sbp(1),sbp(2),[3*sbp(2)+5,3*sbp(2)+6,4*sbp(2)+5,4*sbp(2)+6])
                hold off;
                NpointsPcalc=size(thisPhiNM,1);
                p = hist3(thisPhiNM,{cfg.gridPhi cfg.gridPhi})/NpointsPcalc;
                imagesc(cfg.gridPhi,cfg.gridPhi,p.',[0 2/cfg.Nbins]);
                axis equal
                Cmap=colormap(gray(100));
                Cmap(1:100,:) = Cmap(100:-1:1,:);
                colormap(Cmap);
                colorbar
                set(gca,'YDir','normal')
%                 hold on;
%                 plot(thisPhiNM(:,1),thisPhiNM(:,2),'k.','markersize',0.5);
                xlim([-0.01 1.01 ]*TwoPi)
                ylim([-0.01 1.01 ]*TwoPi)
                axis square
                xlabel([n,'\phi_1'])
                ylabel([m,'\phi_2'])
                set(gca,'xtick',histPhiTick,'ytick',histPhiTick);
                DPformat_ticks_ext_plus(gca,histPhiTickLabelN, histPhiTickLabelM);
                grid on;
                hold off;
                    
                clear thisPhiNM thisDphi;
            end
            
            
            switch method
                
                case 1 %trial
                    
                    
                    figure(h)
                    
                    %Plot this x
                    subplot(sbp(1),sbp(2),sbpPos{1})
                    hold on;
                    for iP=1:2;
                        htemp=findobj(gca,'edgecolor',SignalColors{iP});
                        delete(htemp);
                        hx(iP)=patchline(cfg.time,x(:,iP,iT),'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',1,'linewidth',1);
                    end
                    hold off;
                    
                    %Plot this phi
                    subplot(sbp(1),sbp(2),sbpPos{2})
                    hold on;
                    for iP=1:2;
                        htemp=findobj(gca,'edgecolor',SignalColors{iP});
                        delete(htemp);
                        hphi(iP)=patchline(cfg.timeCut,phiNM(:,iP,iT)/TwoPi,'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',1,'linewidth',1);
                    end
                    hold off;
                    
                    %Plot this Dphi
                    subplot(sbp(1),sbp(2),sbpPos{3})
                    delete(hDphi);
                    hold on;
                    hDphi=patchline(cfg.timeCut,Dphi(:,iT)/pi,'edgecolor','k','facecolor','k','edgealpha',1,'linewidth',1);
                    hold off;
                    
                    %Highlight trial in the measures' plot
                    subplot(sbp(1),sbp(2),sbpPos{4});
                    hold on;
                    iiM=0;
                    xticklabel={};
                    %For each measure we have calculated...
                    if exist('hM','var')
                        delete(hM);
                        clear hM;
                    end
                    for iM = 1:Nmeasures;
                        %...and for each of the sumbmeasures of this measure...
                        for iSubM = 1:NmeasPmeas(measInds(iM));
                            %...plot
                            iiM=iiM+1;
                            hold on;
                            hM(iiM)=plot(iiM,PC.(MeasNames{measInds(iM),iSubM})(iT),'o','markeredgecolor',MeasColors{measInds(iM),iSubM}(1,:),'markersize',10);
                            hold on;
                        end
                    end
                    if exist('htext','var')
%                         htemp=findobj(htext,'String',['\fontsize{14}\bfTrial: ',num2str(iT),'/',num2str(cfg.Ntr)]);
%                         delete(htemp);
                          delete(htext);
                          clear htext;
                    end
                    htext=annotation('textbox',[0.475,0.9,0.1,0.1],'String',['\fontsize{14}\bfTrial: ',num2str(iT),'/',num2str(cfg.Ntr)] ,'color','r','LineStyle','none','FitBoxToText','on');
                    hold off;
                    
                    
                case 2 %trialTime
                    
                    figure(h)
                    
                    %Plot this x
                    subplot(sbp(1),sbp(2),sbpPos{1})
                    hold on;
                    for iP=1:2;
                        htemp=findobj(gca,'edgecolor',SignalColors{iP});
                        delete(htemp);
                        hx(iP)=patchline(cfg.time,x(:,iP,iT),'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',1,'linewidth',1);
                    end
                    hold off;
                    
                    %Plot this phi
                    subplot(sbp(1),sbp(2),sbpPos{2})
                    hold on;
                    for iP=1:2;
                        htemp=findobj(gca,'edgecolor',SignalColors{iP});
                        delete(htemp);
                        hphi(iP)=patchline(cfg.timeCut,phiNM(:,iP,iT)/TwoPi,'edgecolor',SignalColors{iP},'facecolor',SignalColors{iP},'edgealpha',1,'linewidth',1);
                    end
                    hold off;
                    
                    %Plot this Dphi
                    subplot(sbp(1),sbp(2),sbpPos{3})
                    delete(hDphi);
                    hold on;
                    hDphi=patchline(cfg.timeCut,Dphi(:,iT)/pi,'edgecolor','k','facecolor','k','edgealpha',1,'linewidth',1);
                    hold off;
                    
                    %Plot this trials measures' results
                    subplot(sbp(1),sbp(2),sbpPos{4});
                    hold on
                    iiM=0;
                    if exist('hM','var')
                        delete(hM);
                        clear hM;
                    end
                    measLegend={};
                    %For each measure we have calculated...
                    for iM = 1:Nmeasures;
                        %...and for each of the sumbmeasures of this measure...
                        for iSubM = 1:NmeasPmeas(measInds(iM));
                            %...plot
                            iiM=iiM+1;
                            hM(:,iiM)=plot(cfg.timeOut,PC.(MeasNames{measInds(iM),iSubM})(:,iT),marker,'color',MeasColors{measInds(iM),iSubM}(1,:),'markersize',10,'linewidth',2);
                            hold on;
                            measLegend=[measLegend;MeasNames{measInds(iM),iSubM}];
                        end
                    end
                    if exist('htext','var')
%                         htemp=findobj(htext,'String',['\fontsize{14}\bfTrial: ',num2str(iT),'/',num2str(cfg.Ntr)]);
%                         delete(htemp);
                        delete(htext)
                        clear htext
                    end
                    htext=annotation('textbox',[0.5,0.9,0.1,0.1],'String',['\fontsize{14}\bfTrial: ',num2str(iT),'/',num2str(cfg.Ntr)] ,'color','r','LineStyle','none','FitBoxToText','on');
                    hold off;
                    
                    
                otherwise
                    
                    if (cfg.Nwin==1)
                        
                        figure(h)
                        
                        subplot(sbp(1),sbp(2),sbpPos{1})
                        hold on;
                        if exist('hlx','var')
                            delete(hlx);
                            clear hlx;
                        end
                        hlx=plot([cfg.timeCalc(iT) cfg.timeCalc(iT)],[xmin xmax],'r','linewidth',0.5);
                        hold off;
                        
                        subplot(sbp(1),sbp(2),sbpPos{2})
                        hold on;
                        if exist('hlphi','var')
                            delete(hlphi);
                            clear hlphi;
                        end
                        hlphi=plot([cfg.timeCalc(iT) cfg.timeCalc(iT)],[-1 1],'r','linewidth',0.5);
                        hold off;
                        
                        subplot(sbp(1),sbp(2),sbpPos{3})
                        hold on;
                        if exist('hlDphi','var')
                            delete(hlDphi);
                            clear hlDphi;
                        end
                        hlDphi=plot([cfg.timeCalc(iT) cfg.timeCalc(iT)],[-1 1],'r','linewidth',0.5);
                        hold off;
                        
                        subplot(sbp(1),sbp(2),sbpPos{4})
                        hold on;
                        if exist('hlM','var')
                            delete(hlM);
                            clear hlM;
                        end
                        hlM=plot([cfg.timeCalc(iT) cfg.timeCalc(iT)],[yMlim(1) yMlim(2)],'r','linewidth',0.5);
                        xlabel(['\fontsize{14}\bfTime: ',num2str(cfg.timeCut(iT)), ' sec'],'color','r');
                        hold off;
                        
                    else
                        
                        sx = [cfg.timeWinEdgs(iT,1) cfg.timeWinEdgs(iT,2) cfg.timeWinEdgs(iT,2) cfg.timeWinEdgs(iT,1)];
                        sxM = [cfg.timeCalc(iT)-cfg.Ts cfg.timeCalc(iT)+cfg.Ts cfg.timeCalc(iT)+cfg.Ts cfg.timeCalc(iT)-cfg.Ts];
                        
                        figure(h)
                        
                        subplot(sbp(1),sbp(2),sbpPos{1})
                        hold on;
                        if exist('hsx','var')
                            delete(hsx);
                            clear hsx;
                        end
                        hsx=patch(sx,syx,C,'r','edgecolor','r','facealpha',0.1,'edgealpha',1);
                        hold off;
                        
                        
                        subplot(sbp(1),sbp(2),sbpPos{2})
                        hold on;
                        if exist('hsphi','var')
                            delete(hsphi);
                            clear hsphi;
                        end
                        hsphi=patch(sx,sy,C,'r','edgecolor','r','facealpha',0.1,'edgealpha',1);
                        hold off;
                        
                        subplot(sbp(1),sbp(2),sbpPos{3})
                        hold on;
                        if exist('hsDphi','var')
                            delete(hsDphi);
                            clear hsDphi;
                        end
                        hsDphi=patch(sx,sy,C,'r','edgecolor','r','facealpha',0.1,'edgealpha',1);
                        hold off;
                        
                        subplot(sbp(1),sbp(2),sbpPos{4})
                        hold on;
                        if exist('hsM','var')
                            delete(hsM);
                            clear hsM;
                        end
                        hsM=patch(sxM,syM,C,'r','edgecolor','r','facealpha',0.1,'edgealpha',0.1);
                        xlabel(['\fontsize{14}\bfTime: ',num2str(cfg.timeWinEdgs(iT,1)),' - ',num2str(cfg.timeWinEdgs(iT,2)), ' sec'],'color','r');
                        hold off;
                    end
            end %switch statement: plot each frame according to method
            

            % Get frame as an image
            f = getframe(gcf);
            
            
            % Create a colormap for the first frame. For the rest of the frames,
            % use the same colormap
            if iF == 1
                [height, width] = size(rgb2ind(f.cdata, 256, 'nodither'));
                mov = zeros(height, width, 1, Nframes, 'uint8');
                [mov(:,:,1,iF), map] = rgb2ind(f.cdata, 256, 'nodither');
            else
                mov(:,:,1,iF) = rgb2ind(f.cdata, map, 'nodither');
            end
            
        end%for loop: generate each frame
        
        
        % Create animated GIF
        save([filename,'movie','_',n,':',m,'.mat'],'mov','map')
        imwrite(mov, map, [filename,'_',n,':',m,'.gif'], 'DelayTime', 1, 'LoopCount', inf);
        
        
    end %if statement: make animation, separatelly for method==2
    
    hout(iNM)=h;
    
    
    %...if there is statistics, plot it too...
    statsRes=statsResCell{nm(1),nm(2)};
    if ~isempty(statsRes)
        
        
        
        switch method
            
            %case 1
                
            case 2
                
                hstats=figure('Position',[fx fy fwidth fheight],'color','w');
                if cfg.Nout(1)<=50
                    Mmarker = 'ro-';
                    Smarker = 'ko-';
                else
                    Mmarker = 'r';
                    Smarker = 'k';
                end
                iiM=0;
                %For each measure we have calculated...
                for iM = 1:Nmeasures;
                    %...and for each of the sumbmeasures of this measure...
                    for iSubM = 1:NmeasPmeas(measInds(iM));
                        
                        %...plot
                        iiM=iiM+1;
                        
                        subplot(NsubMeasures,2,2*(iiM-1)+1)
                        hold on;
                        plot(cfg.timeOut,statsRes.(MeasNames{measInds(iM),iSubM}).meanPointStat.mVal,Mmarker,'MarkerFaceColor','r');
                        plot(cfg.timeOut,statsRes.(MeasNames{measInds(iM),iSubM}).meanPointStat.surrVal,Smarker,'MarkerFaceColor','k');
                        ylabel(MeasNames{measInds(iM),iSubM});
                        if (iiM==1)
                            legend({'Measurement','Surrogates'})
                            %legend('boxoff')
                        end
                        if (iiM==NsubMeasures)
                            xlabel('Time (sec)')
                        end
                        axis([cfg.time(1)-cfg.Ts cfg.time(end)+cfg.Ts, -0.01, 1.01]);
                        grid on;hold off;
                        
                        subplot(NsubMeasures,2,2*(iiM-1)+2)
                        hold on;
                        plot(cfg.timeOut,cfg.stats.alpha/cfg.stats.tail*ones(size(cfg.timeOut)),'r');
                        plot(cfg.timeOut,statsRes.(MeasNames{measInds(iM),iSubM}).pMean,Smarker,'MarkerFaceColor','k');
                        if (iiM==1)
                            legend({'\alpha-value','p-value'})
                            %legend('boxoff')
                        end
                        if (iiM==NsubMeasures)
                            xlabel('Time (sec)')
                        end
                        axis([cfg.time(1)-cfg.Ts cfg.time(end)+cfg.Ts, -0.01, 1.01]);
                        grid on;hold off
                        
                    end
                end
                %Save figure
                saveas(hstats,[filename,'_',n,':',m,'_Stats.fig']);
                saveas(hstats,[filename,'_',n,':',m,'_Stats.png']);
                
                
            otherwise 
                
                if Nwin==1
                    
                    hstats=figure('Position',[fx fy fwidth fheight],'color','w');
                    iiM=0;
                    %For each measure we have calculated...
                    for iM = 1:Nmeasures;
                        %...and for each of the sumbmeasures of this measure...
                        for iSubM = 1:NmeasPmeas(measInds(iM));
                            
                            %...plot
                            iiM=iiM+1;
                            
                            subplot(NsubMeasures,2,2*(iiM-1)+1)
                            hold on;
                            plot(cfg.time,statsRes.(MeasNames{measInds(iM),iSubM}).pointStat.mVal,'r');
                            plot(cfg.time,statsRes.(MeasNames{measInds(iM),iSubM}).pointStat.surrVal,'k');
                            ylabel(MeasNames{measInds(iM),iSubM});
                            if (iiM==1)
                                legend({'Measurement','Surrogates'})
                                %legend('boxoff')
                            end
                            if (iiM==NsubMeasures)
                                xlabel('Time (sec)')
                            end
                            axis([cfg.time(1)-cfg.Ts cfg.time(end)+cfg.Ts, -0.01, 1.01]);
                            grid on;hold off;
                            
                            subplot(NsubMeasures,2,2*(iiM-1)+2)
                            plot(cfg.time,cfg.stats.alpha/cfg.stats.tail*ones(size(cfg.time)),'r');
                            hold on;
                            plot(cfg.time,statsRes.(MeasNames{measInds(iM),iSubM}).p,'k');
                            if (iiM==1)
                                legend({'\alpha-value','p-value'})
                                %legend('boxoff')
                            end
                            if (iiM==NsubMeasures)
                                xlabel('Time (sec)')
                            end
                            axis([cfg.time(1)-cfg.Ts cfg.time(end)+cfg.Ts, -0.01, 1.01]);
                            grid on;hold off
                            
                        end
                    end
                    %Save figure
                    saveas(hstats,[filename,'_',n,':',m,'_Stats.fig']);
                    saveas(hstats,[filename,'_',n,':',m,'_Stats.png']);
                    
                    
                else
                    
                    hstats=figure('Position',[fx fy fwidth fheight],'color','w');
                    if cfg.Nout(1)<=50
                        Mmarker = 'ro-';
                        Smarker = 'ko-';
                    else
                        Mmarker = 'r';
                        Smarker = 'k';
                    end
                    iiM=0;
                    %For each measure we have calculated...
                    for iM = 1:Nmeasures;
                        %...and for each of the sumbmeasures of this measure...
                        for iSubM = 1:NmeasPmeas(measInds(iM));
                            
                            %...plot
                            iiM=iiM+1;
                            
                            subplot(NsubMeasures,2,2*(iiM-1)+1)
                            hold on;
                            plot(cfg.timeOut,statsRes.(MeasNames{measInds(iM),iSubM}).pointStat.mVal,Mmarker,'MarkerFaceColor','r');
                            plot(cfg.timeOut,statsRes.(MeasNames{measInds(iM),iSubM}).pointStat.surrVal,Smarker,'MarkerFaceColor','k');
                            ylabel(MeasNames{measInds(iM),iSubM});
                            if (iiM==1)
                                legend({'Measurement','Surrogates'})
                            end
                            if (iiM==NsubMeasures)
                                xlabel('Time (sec)')
                            end
                            axis([cfg.time(1)-cfg.Ts cfg.time(end)+cfg.Ts, -0.01, 1.01]);
                            grid on;hold off;
                            
                            subplot(NsubMeasures,2,2*(iiM-1)+2)
                            hold on;
                            plot(cfg.timeOut,cfg.stats.alpha/cfg.stats.tail*ones(size(cfg.timeOut)),'r');
                            plot(cfg.timeOut,statsRes.(MeasNames{measInds(iM),iSubM}).p,Smarker,'MarkerFaceColor','k');
                            if (iiM==1)
                                legend({'\alpha-value','p-value'})
                            end
                            if (iiM==NsubMeasures)
                                xlabel('Time (sec)')
                            end
                            axis([cfg.time(1)-cfg.Ts cfg.time(end)+cfg.Ts, -0.01, 1.01]);
                            grid on;hold off
                            
                        end%for among sub-measures
                    end %for among measures
                    
                    %Save figure
                    saveas(hstats,[filename,'_',n,':',m,'_Stats.fig']);
                    saveas(hstats,[filename,'_',n,':',m,'_Stats.png']);
                end
                
        end %switch method
        
        
    end %is ~isempty(statsRes)
    
    
    
end %for statement: make figures and optionally animations for each nm calculation
h=hout;



% close all