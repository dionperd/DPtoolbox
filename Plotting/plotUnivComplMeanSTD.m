function [h, hsub]=plotUnivComplMeanSTD(M,std,Measure,Colors, Groups, Conds, Ns,cfg, method)

%M = PLSmatStats.MSE.mean
%std = PLSmatStats.MSE.std
%Measure = 'MSE'
%cfg = cfgMeasure
%method='STDerr'


Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);
Nch = size(M{1},2); 

%Make legend
iX=0;
for iG=1:Ng;
    for iC=1:Nc;
        iX=iX+1;
        leg{iX} = [Groups{iG},'-',Conds{iC}];
    end
end
Ndof = Ng*Nc;


if strcmpi(Measure,'MSE')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.MSE.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'MSE';
    
elseif strcmpi(Measure,'MSEn')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.MSEn.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'MSEn';
    
elseif strcmpi(Measure,'LZ')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.LZ.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'LZ';
    
elseif strcmpi(Measure,'LZn')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.LZn.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'LZn';
    
elseif strcmpi(Measure,'CMSEn')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.CMSEn.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'CMSEn';
    
elseif strcmpi(Measure,'CMSE')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.CMSE.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'CMSE';
    
elseif strcmpi(Measure,'LZ')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.LZ.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'LZ';
    
elseif strcmpi(Measure,'LZn')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.LZn.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = 'LZn';
    
elseif strcmpi(Measure,'STD')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.MSE.scales/cfg.fs;
    y=1:Nch;
    xlbl = 'scales (sec)';
    ylbl = '\sigma';
    
elseif strcmpi(Measure,'logF')
    dataMode = 2;
    iCh = 26; %Cz
    x=log(cfg.DFA.scales/cfg.fs);
    y=1:Nch;
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'log(F)';
    
elseif strcmpi(Measure,'logV')
    dataMode = 2;
    iCh = 26; %Cz
    x=log(cfg.VGR.scales/cfg.fs);
    y=1:Nch;
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'log(VAR)';
    
elseif strcmpi(Measure,'logP')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.SP.logf;
    y=1:Nch;
    xlbl = 'log(f) (log(Hz))';
    ylbl = 'log(P)';
    
elseif strcmpi(Measure,'P')
    dataMode = 2;
    iCh = 26; %Cz
    x=cfg.SP.f(2:end);
    y=1:Nch;
    xlbl = 'f (Hz)';
    ylbl = 'P';
    
elseif strcmpi(Measure,'PS')
    dataMode = 1;
    x=[1:Nch].';
    y=[];
    xlbl = 'Channels';
    ylbl = 'Power slope';
elseif strcmpi(Measure,'VS')
    dataMode = 1;
    x=[1:Nch].';
    y=[];
    xlbl = 'Channels';
    ylbl = 'Variogram slope';
elseif strcmpi(Measure,'DOF')
    dataMode = 1;
    x=[1:Nch].';
    y=[];
    xlbl = 'Channels';
    ylbl = 'DoF';
elseif strcmpi(Measure,'H')
    dataMode = 1;
    x=[1:Nch].';
    y=[];
    xlbl = 'Channels';
    ylbl = 'generalized Hurst exponent';
elseif strcmpi(Measure,'mMSE')
    dataMode = 1;
    x=[1:Nch].';
    y=[];
    xlbl = 'Channels';
    ylbl = 'mean MSE';  
elseif strcmpi(Measure,'mMSEn')
    dataMode = 1;
    x=[1:Nch].';
    y=[];
    xlbl = 'Channels';
    ylbl = 'mean normalized MSE';      
end


if strcmpi(method,'STDerr') 
 %plotting mean / standard error ratios maps
    
    figName = [Measure,', mean / standard error'];
    
    maxVal = -inf;
    minVal = inf;
    %Calculate mean / standard error ratio
    for iG=1:Ng;
        for iC=1:Nc;
            xM{iC,iG} = M{iG,iC} ./ (std{iG,iC}/sqrt(Ns(iG))) ;
            minVal = min(minVal,min(min(xM{iC,iG})));
            maxVal = max(maxVal,max(max(xM{iC,iG})));
        end
    end
    
    if (dataMode==1)
        [h,hsub]=topoplot(xM,minVal,maxVal,leg,Ndof,Ng,Nc,figName);
    else
        [h,hsub]=elecsImage(x,y,xM,minVal,maxVal,xlbl,leg,Ndof,Ng,Nc,figName);
        
    end
    
elseif strcmpi(method,'MeanStdErr') 
    
    figName = [Measure,', mean and standard error intervals'];
    
    minVal = inf;
    maxVal = -inf;
    %Calculate standard error intervals
    if (dataMode==1)
        for iG=1:Ng;
            for iC=1:Nc;
                stdErr = std{iG,iC}/sqrt(Ns(iG));
                stdErr = stdErr(:);
                xM{iC,iG} = M{iG,iC}(:);
                ul{iC,iG} = xM{iC,iG} + stdErr;
                ll{iC,iG} = xM{iC,iG} - stdErr;
                minVal = min(minVal,min(ll{iC,iG}(:)));
                maxVal = max(maxVal,max(ll{iC,iG}(:)));
            end
        end
        [h,hsub]=plotLineConfInterv(x,xM,ll,ul,minVal,maxVal,[],xlbl,ylbl,Colors,Groups,Conds,Ng,Nc,figName,dataMode);
        
    else
        for iG=1:Ng;
            for iC=1:Nc;
                stdErr = std{iG,iC}(:,iCh)/sqrt(Ns(iG));
                xM{iC,iG} = M{iG,iC}(:,iCh);
                ul{iC,iG} = M{iG,iC}(:,iCh) + stdErr;
                ll{iC,iG} = M{iG,iC}(:,iCh) - stdErr;
                minVal = min(minVal,min(ll{iC,iG}(:)));
                maxVal = max(maxVal,max(ll{iC,iG}(:)));
            end
        end
        [h,hsub]=plotLineConfInterv(x,xM,ll,ul,minVal,maxVal,iCh,xlbl,ylbl,Colors,Groups,Conds,Ng,Nc,figName,dataMode);
        
    end
else
    %plotting mean maps
    
    figName = [Measure,', mean / standard error'];
    
    maxVal = -inf;
    minVal = inf;
    %Calculate mean / standard error ratio
    for iG=1:Ng;
        for iC=1:Nc;
            xM{iC,iG} = M{iG,iC} ;
            minVal = min(minVal,min(min(xM{iC,iG})));
            maxVal = max(maxVal,max(max(xM{iC,iG})));
        end
    end
    
    if (dataMode==1)
        [h,hsub]=topoplot(xM,minVal,maxVal,leg,Ndof,Ng,Nc,figName);
    else
        [h,hsub]=elecsImage(x,y,xM,minVal,maxVal,xlbl,leg,Ndof,Ng,Nc,figName);
        
    end
end





function [h,hsub]=elecsImage(x,Ch,M,minVal,maxVal,xlbl,leg,Ndof,Ng,Nc,figName)

%19 channels, Cz=10:
%Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
%ChLab = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

%58 channels, Cz=26:
%Ch=[1:58]; %electrodes
ChLab = { 'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4',...
          'F6', 'F8', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'T7', 'C5', 'C3',...
          'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2',...
          'CP4', 'CP6', 'TP8', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8',...
          'PO7', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO8', 'O1', 'Oz', 'O2'};
      
h=figure('Name',figName,'units','normalized','outerposition',[0 0 1 1]);
for iX=1:Ndof;
    hsub(iX)=subplot(Ng,Nc,iX);
    imagesc(x,Ch,M{iX}.',[minVal,maxVal]);
    set(gca,'YDir','reverse','ytick',Ch,'yticklabel',ChLab);%  
    axis tight;
    xlabel(xlbl);ylabel('Channels');
    title(leg{iX})
    colorbar
    hold off;
end

function [h,hsub]=plotLineConfInterv(x,M,ll,ul,minVal,maxVal,iCh,xlbl,ylbl,Colors,Groups,Conds,Ng,Nc,figName,dataMode)
%[h,hsub]=plotLineConfInterv(x,M,ll,ul,minVal,maxVal,iCh,xlbl,ylbl,Colors,leg,Ndof,Ng,Nc,figName)
%19 channels, Cz=10:
%Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
%ChLab = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

%58 channels, Cz=26:
%Ch=[1:58]; %electrodes
ChLab = { 'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4',...
    'F6', 'F8', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'T7', 'C5', 'C3',...
    'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2',...
    'CP4', 'CP6', 'TP8', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8',...
    'PO7', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO8', 'O1', 'Oz', 'O2'};

h=figure('Name',figName,'units','normalized','outerposition',[0 0 1 1]);
for iC=1:Nc;
    hsub(iC)=subplot(Nc,1,iC);
    for iG=1:Ng;
        iX = Nc*(iG-1)+iC;
        hold on;
        if dataMode==1
            errorbar(x,M{iX},M{iX}-ll{iX},ul{iX}-M{iX},'color',Colors{iX},'linewidth',1,'linestyle','none','marker','.');
        else
            p=patch([x; x(end:-1:1)].',[ul{iX}; ll{iX}(end:-1:1)].',ones(1,2*length(x)),'edgecolor',Colors{iX},'facecolor',Colors{iX},'facealpha',0.25,'edgealpha',0);
            hasbehavior(p,'legend',false);
            plot(x,M{iX},'color',Colors{iX},'LineStyle','-','linewidth',2);
        end
    end
    set(gca,'ylim',[minVal,maxVal]);
    xlabel(xlbl);ylabel(ylbl);
    if isempty(iCh)
        title(Conds{iC})
    else
        title([Conds{iC},', ',ChLab{iCh}])
    end
    legend(Groups)
    hold off;
end
% for iX=1:Ndof;
%     hsub(iX)=subplot(Ng,Nc,iX);
%     p=patch([x; x(end:-1:1)].',[ul{iX}; ll{iX}(end:-1:1)].',ones(1,2*length(x)),'edgecolor',Colors{iX},'facecolor',Colors{iX},'facealpha',0.25,'edgealpha',0);
%     hasbehavior(p,'legend',false);
%     plot(x,M{iX},'color',Colors{iX},'LineStyle','-','linewidth',2);
%     set(gca,'ylim',[minVal,maxVal]);
%     xlabel(xlbl);ylabel(ylbl);
%     if isempty(iCh)
%         title(leg{iX})
%     else
%         title([leg{iX},', ',ChLab{iCh}])
%     end
%     hold off;
% end

function h=plotElectrodes(x,LVm,LVll,LVul,xlbl,ylbl,Colors,LineStyle,sbplt)

%19 channels, Cz=10:
%Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
%ChLab = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

%58 channels, Cz=26:
Ch=[1:58]; %electrodes
ChLab = { 'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4',  'F7',  'F5',  'F3',  'F1',  'Fz',  'F2',  'F4',...
           'F6',  'F8', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6',  'T7',  'C5',  'C3',...
           'C1',  'Cz',  'C2',  'C4',  'C6',  'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2',...
          'CP4', 'CP6', 'TP8',  'P7',  'P5',  'P3',  'P1',  'Pz',  'P2', ' P4',  'P6',  'P8',...
          'PO7', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO8',  'O1',  'Oz', 'O2'};
      
% Nc=length(Ch);
%SubPlots=[2,4,6:20,22,24];

iC=26;

minVal = min(min(LVll));
maxVal = max(max(LVul));
%for iC=1:Nc;
    subplot(sbplt(1),sbplt(2),sbplt(3));
    hold on;
%     p=patch([x; x(end:-1:1)],[LVul(:,Ch(iC)); LVll(end:-1:1,Ch(iC))],ones(1,2*length(x)),'edgecolor',Color,'facecolor',Color,'facealpha',0.25,'edgealpha',0);
%     hasbehavior(p,'legend',false);
%     plot(x,LVm(:,Ch(iC)),'color',Colors{iX},'LineStyle','-','linewidth',2);
    errorbar(x,LVm(:,Ch(iC)),LVm(:,Ch(iC))-LVll(:,Ch(iC)),LVul(:,Ch(iC))-LVm(:,Ch(iC)),'color',Colors{iX},'linewidth',1,'linestyle','none','marker','.')
    h=gca;
    set(h,'ylim',[minVal,maxVal]);
    grid on;
    title(ChLab{iC})
    xlabel(xlbl);ylabel(ylbl);
    
%end
%h=gcf;


function  [h, hsub]=topoplot(M,minVal,maxVal,leg,Ndof,Ng,Nc,figName)

% %Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
% 
% chx = [ 53.0000   53.0000  107.0000  121.0000  125.0000  122.0000  114.0000  204.0000  199.0000...
%         198.0000  200.0000  206.0000  288.0000  276.0000  271.0000  276.0000  288.0000  341.0000  342.0000].';
%     
% chy = [ 148.0000  239.0000   71.0000  133.0000  192.0000  255.0000  313.0000   51.0000  119.0000...
%         191.0000  259.0000  332.0000   75.0000  128.0000  190.0000  253.0000  305.0000  147.0000  231.0000].';
% LV = LV(Ch);
% Ch=[1:19]; %electrodes

%For 58 electrodes:
chx = [ 50.9956,   57.9605,   54.4781,   86.9810,   86.9810,  107.8757,  114.8406,  120.6447,  124.1272,  126.4488, 125.2880,  121.8056,  118.3231,  111.3582,  155.4693,  158.9518,  163.5950,  162.4342,  162.4342,  161.2734, 156.6301,  200.7412,  200.7412,  198.4196,  197.2588,  197.2588,  196.0980,  198.4196,  198.4196,  207.7061, 246.0132,  246.0132,  239.0482,  236.7266,  235.5658,  237.8874,  241.3699,  246.0132,  248.3348,  283.1594, 279.6769,  273.8728,  271.5512,  269.2295,  272.7120,  273.8728,  284.3202, 287.8026,  316.8231,  312.1798, 301.7325,  306.3757,  304.0541, 313.3406,  321.4664,  340.0395,  341.2003,  342.3611].';

chy = [  146.3710  193.3387  237.5968,  135.5323  250.2419   70.5000  100.3065  131.0161  163.5323  190.6290 223.1452  254.7581  283.6613  316.1774   88.5645  123.7903  157.2097  190.6290  224.9516  258.3710 294.5000   50.6290   85.8548  119.2742  153.5968  191.5323  227.6613  261.0806  297.2097  334.2419 52.4355   92.1774  120.1774  152.6935  190.6290  223.1452  255.6613  291.7903  323.4032   75.0161 104.8226  125.5968  158.1129  190.6290  220.4355  252.0484  278.2419  305.3387  103.0161  136.4355 166.2419  187.0161  214.1129  243.9194  272.8226  148.1774  187.9194  229.4677 ]';


h=figure('Name',figName,'units','normalized','outerposition',[0 0 1 1]);
for iX=1:Ndof;
    hsub(iX)=subplot(Ng,Nc,iX);
    %[h, z,map]=eegplot(mag,ch,unorm,ch_disp,method,color_res,minVal,maxVal,ColBar)
    [dummy,dummy1,dummy2]=eegplot(double(M{iX}(:)),[chx, chy],1,1,'nearest',[],minVal,maxVal,1);
    title(leg{iX})
    hold off;
end
